##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## AI usage telemetry — a flat, append-only, frozen-priced event table.
##
## One row per AI call lands in a single common SQLite store
## (etc/ai_telemetry.sqlite, carrying a user_email column) that ships to
## BigOmics servers for downstream analytics. Two ingestion paths feed it:
##
##   * ai_telemetry_collect_chat()  — PULL: scans assistant turns in a user's
##     chat store (sessions-<user>.sqlite) and appends one event per turn.
##     Idempotent via INSERT OR IGNORE on event_id = "<session>#t<turn_idx>".
##   * ai_telemetry_record()        — PUSH: one-shot sink for reports / cards /
##     wgcna summaries that generate AI text outside the chat flow.
##
## Token parsing is delegated wholesale to omicsai / omicsagentovi:
##   omicsagentovi::turn_from_payload() reconstructs an ellmer::Turn from the
##   stored payload, and omicsai::omicsai_turn_usage() returns a typed S7
##   <omicsai_usage> (slots read with @, never $). We never re-implement the
##   provider-specific json.usage path parsing here.

## ---------------------------------------------------------------------------
## Default paths
##
## ETC is the app-wide config dir (components/app/R/global.R), but it is unbound
## when this file is sourced standalone (tests, CLI). Resolve lazily so the
## signatures stay valid at source time and fall back to a relative "etc".

.ai_tel_default_db <- function() {
  file.path(if (exists("ETC")) ETC else "etc", "ai_telemetry.sqlite")
}

.ai_tel_default_pricing <- function() {
  file.path(if (exists("ETC")) ETC else "etc", "model_pricing.yml")
}

## ---------------------------------------------------------------------------
## .2 — DB layer: schema + open

#' Create the ai_usage_events table and set WAL pragmas
#'
#' Idempotent (CREATE TABLE IF NOT EXISTS). busy_timeout is essential here —
#' unlike a per-session chat store this common db has many concurrent writers.
#' @param con An open DBI connection.
#' @return invisible(con)
.ai_tel_ensure_schema <- function(con) {
  DBI::dbExecute(con, "PRAGMA journal_mode = WAL")
  DBI::dbExecute(con, "PRAGMA synchronous = NORMAL")
  DBI::dbExecute(con, "PRAGMA busy_timeout = 5000")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS ai_usage_events (
      event_id        TEXT PRIMARY KEY,
      source          TEXT,
      session_id      TEXT,
      turn_idx        INTEGER,
      user_email      TEXT,
      model           TEXT,
      created_at      REAL,
      recorded_at     REAL,
      tool_count      INTEGER,
      input_fresh     INTEGER,
      input_cached    INTEGER,
      output          INTEGER,
      reasoning       INTEGER,
      total           INTEGER,
      cache_hit_rate  REAL,
      approx_price_usd REAL,
      pricing_version TEXT
    )")
  invisible(con)
}

#' Open the common telemetry db (read-write) and ensure its schema
#'
#' Caller owns the connection and must DBI::dbDisconnect() it.
#' @param db_path Path to the SQLite file. Created if missing.
#' @param create When TRUE, ensure the parent dir + schema exist.
#' @return An open DBI connection.
ai_telemetry_open <- function(db_path = .ai_tel_default_db(), create = TRUE) {
  if (create) {
    dir <- dirname(db_path)
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  if (create) .ai_tel_ensure_schema(con)
  con
}

## ---------------------------------------------------------------------------
## .1 — Pricing table + resolver

# Parsed-YAML cache, keyed by "<path>@<mtime>" so a changed file reloads but
# repeated calls within a session do not re-read from disk.
.ai_tel_pricing_cache <- new.env(parent = emptyenv())

#' Load (and cache) the pricing YAML.
#' @return The parsed list, or NULL when the file is missing / unreadable.
.ai_tel_load_pricing <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  key <- paste0(path, "@", as.numeric(file.info(path)$mtime))
  cached <- .ai_tel_pricing_cache[[key]]
  if (!is.null(cached)) return(cached)
  parsed <- tryCatch(yaml::read_yaml(path), error = function(e) NULL)
  if (!is.null(parsed)) .ai_tel_pricing_cache[[key]] <- parsed
  parsed
}

#' Resolve per-1M-token rates for a model string
#'
#' Glob-matches `model` against the YAML keys. Precedence: an exact key wins,
#' then the most-specific matching glob (longest non-"*" pattern), then the
#' "*" catch-all. Matches both the session string ("openai:gpt-5.4-nano") and
#' the dated per-turn model ("gpt-5.4-nano-2026-03-17").
#' @param model Character model string.
#' @param path  Path to model_pricing.yml.
#' @return list(input, cached, output, version), or NULL when unresolved.
ai_telemetry_pricing <- function(model, path = .ai_tel_default_pricing()) {
  spec <- .ai_tel_load_pricing(path)
  if (is.null(spec) || is.null(spec$models) || is.null(model) || is.na(model)) {
    return(NULL)
  }
  keys <- names(spec$models)
  hit <- NULL

  # 1. exact match
  if (model %in% keys) {
    hit <- model
  } else {
    # 2. most-specific glob (exclude bare "*"; longest pattern wins)
    globs <- setdiff(keys, "*")
    matched <- globs[vapply(globs, function(k) {
      grepl(utils::glob2rx(k), model)
    }, logical(1))]
    if (length(matched)) {
      hit <- matched[which.max(nchar(matched))]
    } else if ("*" %in% keys) {
      # 3. catch-all
      hit <- "*"
    }
  }
  if (is.null(hit)) return(NULL)

  rates <- spec$models[[hit]]
  list(
    input   = rates$input,
    cached  = rates$cached,
    output  = rates$output,
    version = spec$version
  )
}

#' Compute the frozen USD price for one usage record
#'
#' price = input_fresh/1e6*input + input_cached/1e6*cached +
#'         (output + reasoning)/1e6*output
#' NA token counts are treated as 0. Returns NA when pricing is unresolved.
#' @param usage list/omicsai_usage-like with input_fresh, input_cached, output,
#'   reasoning fields.
#' @param model Model string (used to resolve pricing when not supplied).
#' @param pricing Pre-resolved pricing list (default: resolve from model).
#' @return numeric(1) USD, rounded to 6 dp.
ai_telemetry_price <- function(usage, model, pricing = ai_telemetry_pricing(model)) {
  if (is.null(pricing)) return(NA_real_)
  known <- function(x) !(is.null(x) || length(x) != 1L || is.na(x))
  # No token detail at all (e.g. the report future-worker path) -> price is
  # unknown, NOT zero: a $0 row would understate cost in downstream analytics.
  if (!any(known(usage$input_fresh), known(usage$input_cached),
           known(usage$output), known(usage$reasoning))) {
    return(NA_real_)
  }
  z <- function(x) if (is.null(x) || length(x) != 1L || is.na(x)) 0 else as.numeric(x)
  price <-
    z(usage$input_fresh)  / 1e6 * z(pricing$input) +
    z(usage$input_cached) / 1e6 * z(pricing$cached) +
    (z(usage$output) + z(usage$reasoning)) / 1e6 * z(pricing$output)
  round(price, 6)
}

## ---------------------------------------------------------------------------
## .3 — Turn parser

#' Parse a stored assistant turn payload into a flat usage row
#'
#' Reconstructs an ellmer::Turn via omicsagentovi::turn_from_payload() and reads
#' tokens with omicsai::omicsai_turn_usage() (passing the session model so the
#' provider — hence cache/reasoning JSON paths — resolves correctly).
#'
#' Never throws on a malformed turn: on failure it returns NA token fields with
#' a best-effort tool_count so collect_chat() can still record the event.
#'
#' @param payload_json The raw JSON string from turns.payload_json, OR an
#'   already-parsed payload list.
#' @param session_model The model from the sessions row (used for pricing and
#'   provider resolution).
#' @return list(model, input_fresh, input_cached, output, reasoning, total,
#'   cache_hit_rate, tool_count)
.ai_tel_parse_turn <- function(payload_json, session_model = NULL) {
  na_row <- function(model = session_model, tool_count = 0L) {
    list(
      model          = model,
      input_fresh    = NA_integer_,
      input_cached   = NA_integer_,
      output         = NA_integer_,
      reasoning      = NA_integer_,
      total          = NA_integer_,
      cache_hit_rate = NA_real_,
      tool_count     = tool_count
    )
  }

  payload <- if (is.character(payload_json)) {
    tryCatch(jsonlite::fromJSON(payload_json, simplifyVector = FALSE),
             error = function(e) NULL)
  } else {
    payload_json
  }
  if (!is.list(payload)) return(na_row())

  # tool_count is cheap + robust; compute it before the riskier usage path.
  tool_count <- tryCatch(
    sum(vapply(payload$contents, function(c) identical(c$type, "tool_request"),
               logical(1))),
    error = function(e) 0L
  )
  tool_count <- as.integer(tool_count %||% 0L)

  model <- payload$json$model %||% session_model

  # JSON storage loses the numeric-NA distinction on `tokens` (stored as "NA"
  # strings). Mirror session-store.R restore: flatten + coerce back to numeric
  # NA before turn_from_payload(). See omicsagentovi/R/session-store.R ~L1116.
  tok <- payload$tokens
  if (!is.null(tok)) {
    flat <- suppressWarnings(as.numeric(unlist(tok)))
    payload$tokens <- if (length(flat)) flat else NULL
  }

  u <- tryCatch({
    turn <- omicsagentovi::turn_from_payload(payload)
    omicsai::omicsai_turn_usage(turn, model = session_model)
  }, error = function(e) NULL)

  if (is.null(u)) return(na_row(model = model, tool_count = tool_count))

  list(
    model          = model,
    input_fresh    = u@input_tokens_fresh,
    input_cached   = u@input_tokens_cached,
    output         = u@output_tokens,
    reasoning      = u@reasoning_tokens,
    total          = u@total_tokens,
    cache_hit_rate = u@cache_hit_rate,
    tool_count     = tool_count
  )
}

## ---------------------------------------------------------------------------
## Usage normalisation (shared by record() and the insert helper)

# Coerce an omicsai_usage S7 object / legacy list / NULL into the four token
# fields the pricing + schema use. S7 is accessed with @ (the $ accessor
# silently returns NULL on S7), legacy lists with $, NULL maps to all-NA.
.ai_tel_usage_fields <- function(usage) {
  if (is.null(usage)) {
    return(list(input_fresh = NA_integer_, input_cached = NA_integer_,
                output = NA_integer_, reasoning = NA_integer_,
                total = NA_integer_, cache_hit_rate = NA_real_))
  }
  if (S7::S7_inherits(usage, omicsai::omicsai_usage)) {
    return(list(
      input_fresh    = usage@input_tokens_fresh,
      input_cached   = usage@input_tokens_cached,
      output         = usage@output_tokens,
      reasoning      = usage@reasoning_tokens,
      total          = usage@total_tokens,
      cache_hit_rate = usage@cache_hit_rate
    ))
  }
  # legacy list shape from .omicsai_extract_usage()
  list(
    input_fresh    = usage$input_tokens_fresh,
    input_cached   = usage$input_tokens_cached,
    output         = usage$output_tokens,
    reasoning      = usage$reasoning_tokens,
    total          = usage$total_tokens,
    cache_hit_rate = usage$cache_hit_rate
  )
}

# Single INSERT OR IGNORE for one event row. Returns the number of rows
# actually written (0 when the event_id already exists). Assumes `con` already
# has the schema. NA-safe via as.integer/as.numeric on each field.
.ai_tel_insert_event <- function(con, fields) {
  i <- function(x) if (is.null(x) || length(x) != 1L) NA_integer_ else suppressWarnings(as.integer(x))
  n <- function(x) if (is.null(x) || length(x) != 1L) NA_real_ else suppressWarnings(as.numeric(x))
  s <- function(x) if (is.null(x) || length(x) != 1L || is.na(x)) NA_character_ else as.character(x)

  DBI::dbExecute(con, "
    INSERT OR IGNORE INTO ai_usage_events
      (event_id, source, session_id, turn_idx, user_email, model,
       created_at, recorded_at, tool_count,
       input_fresh, input_cached, output, reasoning, total,
       cache_hit_rate, approx_price_usd, pricing_version)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
    params = list(
      s(fields$event_id), s(fields$source), s(fields$session_id),
      i(fields$turn_idx), s(fields$user_email), s(fields$model),
      n(fields$created_at), n(fields$recorded_at), i(fields$tool_count),
      i(fields$input_fresh), i(fields$input_cached), i(fields$output),
      i(fields$reasoning), i(fields$total),
      n(fields$cache_hit_rate), n(fields$approx_price_usd),
      s(fields$pricing_version)
    ))
  changed <- DBI::dbGetQuery(con, "SELECT changes() AS c")$c
  as.integer(changed %||% 0L)
}

## ---------------------------------------------------------------------------
## .7 — record(): one-shot push sink

#' Record a single one-shot AI usage event (reports / cards / wgcna summaries)
#'
#' Opens the common db, ensures the schema, prices the usage, and writes one row
#' with INSERT OR IGNORE. The connection is closed before returning.
#'
#' @param source    Source tag, e.g. "report", "card", "wgcna_summary".
#' @param session_id Session id (one-shot calls supply their own).
#' @param user_email From auth$email; may be NA.
#' @param model     Model string (used for pricing + provider resolution).
#' @param usage     omicsai_usage S7 / legacy list / NULL.
#' @param turn_idx  Turn index, or NA for one-shot (drives event_id form).
#' @param tool_count Number of tool calls in the turn.
#' @param created_at Call time (POSIXct or numeric).
#' @param db_path    Telemetry db path.
#' @return invisible(event_id)
ai_telemetry_record <- function(source, session_id, user_email, model, usage,
                                turn_idx = NA, tool_count = 0L,
                                created_at = Sys.time(),
                                db_path = .ai_tel_default_db()) {
  # One-shot push events (cards/reports/summaries) have no re-scan, so each
  # record() call is a DISTINCT event. Keying on session_id alone would collapse
  # every generation in a Shiny session into one row (INSERT OR IGNORE drops the
  # rest) — so include the call timestamp (µs) for per-call uniqueness. The chat
  # PULL path stays deterministic ("<session>#t<turn>") to preserve idempotency.
  event_id <- if (is.na(turn_idx)) {
    paste0(source, ":", session_id, ":", sprintf("%.6f", as.numeric(created_at)))
  } else {
    paste0(session_id, "#t", turn_idx)
  }

  uf <- .ai_tel_usage_fields(usage)
  pricing <- ai_telemetry_pricing(model)
  price <- ai_telemetry_price(uf, model, pricing)

  fields <- c(
    list(
      event_id        = event_id,
      source          = source,
      session_id      = session_id,
      turn_idx        = turn_idx,
      user_email      = user_email,
      model           = model,
      created_at      = as.numeric(created_at),
      recorded_at     = as.numeric(Sys.time()),
      tool_count      = tool_count,
      approx_price_usd = price,
      pricing_version = pricing$version
    ),
    uf
  )

  con <- ai_telemetry_open(db_path, create = TRUE)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  .ai_tel_insert_event(con, fields)

  invisible(event_id)
}

## ---------------------------------------------------------------------------
## .5 — collect_chat(): idempotent pull from a user's chat store

#' Ingest a user's assistant chat turns into the common telemetry db
#'
#' Reads sessions-<user>.sqlite read-only (WAL: readers never block writers),
#' joins each assistant turn to its session model, prices it, and appends one
#' event per turn with INSERT OR IGNORE (event_id = "<session>#t<turn_idx>").
#' Fully synchronous and re-runnable — running twice yields the same row count.
#'
#' @param session_dir The user's chats/ directory.
#' @param user_email  From auth$email.
#' @param db_path     Telemetry db path.
#' @return integer count of newly inserted rows.
ai_telemetry_collect_chat <- function(session_dir, user_email,
                                      db_path = .ai_tel_default_db()) {
  chat_db <- .detect_legacy_sessions(session_dir, user_email)
  if (!file.exists(chat_db)) return(0L)

  ro <- DBI::dbConnect(RSQLite::SQLite(), chat_db, flags = RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(ro), add = TRUE)

  rows <- tryCatch(
    DBI::dbGetQuery(ro, "
      SELECT t.session_id, t.turn_idx, t.payload_json, s.model
        FROM turns t
        LEFT JOIN sessions s ON s.session_id = t.session_id
       WHERE t.role = 'assistant'"),
    error = function(e) NULL
  )
  if (is.null(rows) || nrow(rows) == 0L) return(0L)

  con <- ai_telemetry_open(db_path, create = TRUE)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  n_new <- 0L
  for (k in seq_len(nrow(rows))) {
    session_id <- rows$session_id[[k]]
    turn_idx   <- rows$turn_idx[[k]]
    model_in   <- rows$model[[k]]
    parsed <- .ai_tel_parse_turn(rows$payload_json[[k]], session_model = model_in)

    pricing <- ai_telemetry_pricing(parsed$model)
    price <- ai_telemetry_price(parsed, parsed$model, pricing)

    fields <- list(
      event_id         = paste0(session_id, "#t", turn_idx),
      source           = "chat",
      session_id       = session_id,
      turn_idx         = turn_idx,
      user_email       = user_email,
      model            = parsed$model,
      created_at       = NA_real_,
      recorded_at      = as.numeric(Sys.time()),
      tool_count       = parsed$tool_count,
      input_fresh      = parsed$input_fresh,
      input_cached     = parsed$input_cached,
      output           = parsed$output,
      reasoning        = parsed$reasoning,
      total            = parsed$total,
      cache_hit_rate   = parsed$cache_hit_rate,
      approx_price_usd = price,
      pricing_version  = pricing$version
    )
    n_new <- n_new + .ai_tel_insert_event(con, fields)
  }
  as.integer(n_new)
}

## ---------------------------------------------------------------------------
## App-scope helpers: chat-db naming + one-shot legacy migration.
## Single source of truth for the per-user chat filename; the Shiny worker
## (CopilotBoardServer / SessionStore wiring) calls these.

#' Per-user chat-store filename
#'
#' user_dir is already pgx_dir/<email> so the email is filesystem-safe, but
#' defensively neutralise path separators and control characters.
#' @param email User email.
#' @return "sessions-<fs-safe email>.sqlite"
.ai_tel_sessions_db_name <- function(email) {
  if (is.null(email) || length(email) != 1L || is.na(email) || !nzchar(email)) {
    return("sessions.sqlite")
  }
  safe <- gsub("[/\\\\[:cntrl:]]", "_", email)
  paste0("sessions-", safe, ".sqlite")
}

#' Resolve the per-user chat db, migrating a legacy sessions.sqlite once
#'
# DELETE AFTER ~2026-07-23 once dev team migrated (epic omicsplayground-iyp.4).
# AI usage has been internal-only so a single rename pass is safe: if the legacy
# shared sessions.sqlite exists and the per-user target does not, rename it
# (plus its -wal/-shm/-journal sidecars) into place.
#'
#' @param chat_dir The user's chats/ directory.
#' @param email    User email.
#' @return Path to the resolved per-user db.
.detect_legacy_sessions <- function(chat_dir, email) {
  target <- file.path(chat_dir, .ai_tel_sessions_db_name(email))
  legacy <- file.path(chat_dir, "sessions.sqlite")
  if (file.exists(legacy) && !file.exists(target)) {
    file.rename(legacy, target)
    for (ext in c("-wal", "-shm", "-journal")) {
      side_old <- paste0(legacy, ext)
      if (file.exists(side_old)) file.rename(side_old, paste0(target, ext))
    }
  }
  target
}

## ---------------------------------------------------------------------------
## .8 — record_reports(): reconcile helper for AI report slots

#' Record telemetry for all AI report slots that carry usage data
#'
#' Scans \code{pgx$ai} and calls \code{ai_telemetry_record()} once per slot
#' that has a non-NULL \code{usage} field. Uses the slot's \code{created_at}
#' epoch and a dataset-keyed \code{session_id} so the resulting
#' \code{event_id = "report:<dataset_id>:<module>:<created_at>"} is stable
#' across sessions and reloads — INSERT OR IGNORE deduplicates replays.
#'
#' Safe to call multiple times on the same pgx; tolerate \code{user_email = NA}.
#' Errors are logged but never propagate (fire-and-forget).
#'
#' @param pgx       PGX object as a plain list (not reactive).
#' @param user_email User email from auth\$email; may be NA.
#' @param db_path   Telemetry db path.
#' @return Invisible character vector of event ids recorded (may be empty).
ai_telemetry_record_reports <- function(pgx, user_email,
                                        db_path = .ai_tel_default_db()) {
  if (is.null(pgx) || !is.list(pgx)) return(invisible(character(0)))
  ai <- if (is.list(pgx$ai)) pgx$ai else NULL
  if (is.null(ai)) return(invisible(character(0)))

  dataset_id <- if (!is.null(pgx$name) && nzchar(pgx$name %||% "")) pgx$name else ""
  slots      <- names(ai)
  recorded   <- character(0)

  for (module in slots) {
    slot <- ai[[module]]
    if (!is.list(slot)) next
    usage <- slot$usage
    if (is.null(usage)) next

    model      <- usage$model %||% NA_character_
    created_at <- slot$created_at %||% as.numeric(Sys.time())

    # Build a stable, dataset-keyed session_id so the sink's formula produces:
    #   event_id = "report:<dataset_id>:<module>:<created_at_us>"
    session_id <- paste0(dataset_id, ":", module)

    ev_id <- tryCatch(
      ai_telemetry_record(
        source     = "report",
        session_id = session_id,
        user_email = user_email,
        model      = model,
        usage      = usage,
        created_at = structure(created_at, class = c("POSIXct", "POSIXt")),
        db_path    = db_path
      ),
      error = function(e) {
        warning("[ai_telemetry_record_reports] failed for module '", module,
                "': ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (!is.null(ev_id)) recorded <- c(recorded, ev_id)
  }

  invisible(recorded)
}
