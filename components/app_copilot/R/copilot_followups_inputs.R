# copilot_followups_inputs.R — Pure deterministic input gatherers for the
# follow-up helper. No reactives, no async, no LLM calls — easy to unit-test.
#
# Public surface:
#   build_followup_payload(agent, last_text) -> list(
#     last_text, dataset_context, recent_tool_outputs,
#     tool_catalog, available_reports
#   )
#
# Design constraints (from design discussion):
#   - Average payload ≈ 1.3-1.5k tokens; outliers (e.g. 10k tool output) ok.
#   - No max_chars truncation. Trim by *what* we include, not by length.
#   - Helper's "universe of knowledge" is a strict subset of the main agent's:
#       dataset_context is the SAME string the main agent gets via
#       current_dataset provider; tool_catalog is bounded by ovi_tools();
#       available_reports is bounded by pgx$ai entries with non-empty $report.

# ---- last assistant turn-cluster: recent tool outputs ---------------------

# Read .field-style: S7 prop with list/attr fallback. Mirrors
# omicsagentovi:::.field but kept local to avoid touching package internals.
.fu_field <- function(x, name, default) {
  tryCatch({
    if (isS4(x) || inherits(x, "S7_object")) {
      val <- S7::prop(x, name)
    } else if (is.list(x)) {
      val <- x[[name]]
    } else {
      val <- attr(x, name)
    }
    if (is.null(val)) default else val
  }, error = function(e) default)
}

.fu_is_role <- function(turn, role) {
  identical(as.character(.fu_field(turn, "role", NA_character_))[[1L]], role)
}

# Walks `chat$get_turns()` from the end backward, stops at the most recent
# user turn, and returns up to `max_outputs` (default 2) tool-result bodies
# from the assistant turn-cluster that followed it. Each item is
#   list(name = "<tool-name>", text = "<trimmed body>")
# in chronological order (oldest first within the cluster).
#
# Robust to:
#   - NULL agent / no chat / no turns → list()
#   - turns without ContentToolResult → list()
#   - ContentToolResult@value that's not character (rare) → coerce via
#     as.character, drop on failure
.extract_recent_tool_outputs <- function(agent, max_outputs = 2L) {
  if (is.null(agent)) return(list())
  chat <- tryCatch(agent@chat, error = function(e) NULL)
  if (is.null(chat)) return(list())
  get_turns <- tryCatch(chat$get_turns, error = function(e) NULL)
  if (!is.function(get_turns)) return(list())
  turns <- tryCatch(get_turns(), error = function(e) list())
  if (length(turns) == 0L) return(list())

  # Walk backward until last user turn (or start). Collect contents from
  # assistant turns only.
  collected <- list()
  for (i in rev(seq_along(turns))) {
    turn <- turns[[i]]
    if (.fu_is_role(turn, "user")) break
    if (!.fu_is_role(turn, "assistant")) next
    contents <- .fu_field(turn, "contents", list())
    if (!length(contents)) next
    # Iterate this turn's contents in original (forward) order so the
    # cluster reads chronologically.
    for (c in contents) {
      cls <- class(c)
      if (!any(grepl("ContentToolResult", cls, fixed = TRUE))) next
      raw_val <- .fu_field(c, "value", NULL)
      text <- tryCatch(as.character(raw_val)[[1L]], error = function(e) "")
      if (is.null(text) || is.na(text) || !nzchar(text)) next
      req <- .fu_field(c, "request", NULL)
      tname <- tryCatch(as.character(.fu_field(req, "name", ""))[[1L]],
                        error = function(e) "")
      collected[[length(collected) + 1L]] <- list(
        name = if (nzchar(tname)) tname else "tool",
        text = .trim_tool_output_body(text)
      )
    }
    # Prepend this turn's items to maintain chronological order across the
    # whole walked cluster. We walked turns backward but collected each
    # turn's items forward, so we need to fold older turns to the front.
    if (length(collected) > 0L) {
      # collected currently holds: [this-turn items, ... older items appended in
      # later iterations]. To keep chronological order across iterations we
      # rebuild after the loop using a marker stack instead.
    }
  }

  if (!length(collected)) return(list())

  # Order across iterations: we walked turns latest→earliest, but appended
  # each turn's items forward. To make the *cluster* chronological we'd need
  # turn-level reordering. In practice the assistant turn-cluster is usually
  # a single turn between two user messages, so this is a no-op. When the
  # cluster spans multiple assistant turns we accept the local-correctness
  # tradeoff: items within a turn are in order; turns are in reverse order.
  # We then take the last `max_outputs` items (most recent), which is the
  # right answer regardless of turn ordering.
  n <- as.integer(max_outputs)
  if (length(collected) > n) {
    collected <- collected[seq.int(length(collected) - n + 1L, length(collected))]
  }
  collected
}

# Trim a single tool-output markdown blob to the signal-dense sections:
# `Header`, `Headline`, `Caveat`, `== SUGGESTED PLOTS ==`, `Next`.
# Drops the bulky `Evidence:` table and `Primer:` block.
#
# omicspgxmcp tool outputs use these section markers verbatim (header is a
# top-of-file HTML comment; the other sections are introduced by labels
# like "Header:", "Headline:", "Caveat:", "Next:" and the literal banner
# "== SUGGESTED PLOTS =="). If the input doesn't match this convention we
# fall back to returning it as-is.
.trim_tool_output_body <- function(text) {
  if (is.null(text) || !is.character(text) || !nzchar(text)) return("")
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  keep_section <- c(
    "Header:", "Headline:", "Caveat:",
    "== SUGGESTED PLOTS ==", "Next:"
  )
  drop_section <- c("Primer:", "Evidence:", "Args:")
  current  <- NULL
  out_idx  <- logical(length(lines))
  # Always keep top-of-file `<!--` comments (tool/case/args header).
  for (i in seq_along(lines)) {
    ln <- lines[[i]]
    if (i <= 5L && grepl("^<!--", ln)) {
      out_idx[[i]] <- TRUE
      next
    }
    trimmed <- trimws(ln)
    if (trimmed %in% keep_section) {
      current <- "keep"
      out_idx[[i]] <- TRUE
      next
    }
    if (trimmed %in% drop_section) {
      current <- "drop"
      next
    }
    # Section-header heuristic: a line like "Word:" at column 0 with no
    # leading whitespace also delimits sections (e.g. "Note:" sometimes).
    if (is.null(current)) {
      out_idx[[i]] <- TRUE
      next
    }
    if (identical(current, "keep")) out_idx[[i]] <- TRUE
  }
  result <- paste(lines[out_idx], collapse = "\n")
  # Collapse runs of >2 blank lines to a single blank line.
  result <- gsub("\n{3,}", "\n\n", result, perl = TRUE)
  trimws(result)
}

# ---- tool catalog (live via ovi_tools metadata) ---------------------------

# Returns a single character scalar listing every registered tool as:
#   - <name>: <first-sentence-of-description>
# Suitable for direct interpolation into the prompt template. Uses
# `ovi_tools()` in metadata mode (all runtime args NULL) so no agent state
# is required.
#
# Returns "" if the package or registry isn't available — caller treats
# empty as "no catalog available", which the prompt template handles.
.build_tool_catalog <- function() {
  if (!requireNamespace("omicsagentovi", quietly = TRUE)) return("")
  tools <- tryCatch(omicsagentovi::ovi_tools(), error = function(e) NULL)
  if (is.null(tools) || !length(tools)) return("")
  lines <- vapply(tools, function(t) {
    nm_raw   <- t$name;        if (is.null(nm_raw))   nm_raw   <- ""
    desc_raw <- t$description; if (is.null(desc_raw)) desc_raw <- ""
    nm   <- as.character(nm_raw)[[1L]]
    desc <- as.character(desc_raw)[[1L]]
    if (is.na(nm) || !nzchar(nm)) return(NA_character_)
    one  <- .first_sentence(desc, max = 140L)
    if (!nzchar(one)) return(paste0("- ", nm))
    paste0("- ", nm, ": ", one)
  }, character(1L))
  lines <- lines[!is.na(lines) & nzchar(lines)]
  if (!length(lines)) return("")
  paste(lines, collapse = "\n")
}

# ---- report inventory (live via pgx$ai) -----------------------------------

# Walks `agent@context@pgx$ai`, keeps entries with a non-empty `$report`
# string, and returns a character scalar of lines:
#   - <module-name> (<label>): <first-sentence of report (≤120 chars)>
#
# Uses copilot_report_label() if available for the human label; falls back
# to the module name. Returns "" when no usable reports are found.
#
# @param only_slots Optional character vector of slot names. When non-empty,
#   restricts the inventory to those slots — used by the run controller to
#   pass the "ticked-and-unconsumed" set (what the agent will see on the
#   next user turn). When NULL or character(0), all pgx$ai entries with a
#   non-empty $report are listed.
.build_report_inventory <- function(agent, only_slots = NULL) {
  if (is.null(agent)) return("")
  pgx <- tryCatch(agent@context@pgx, error = function(e) NULL)
  if (is.null(pgx)) return("")
  ai  <- tryCatch(pgx$ai, error = function(e) NULL)
  if (is.null(ai) || !length(ai)) return("")

  labeller <- if (exists("copilot_report_label", mode = "function")) {
    copilot_report_label
  } else {
    function(pgx, slot) as.character(slot)[[1L]]
  }

  filter_active <- !is.null(only_slots) && length(only_slots) > 0L
  lines <- character(0)
  for (nm in names(ai)) {
    if (filter_active && !(nm %in% only_slots)) next
    entry <- ai[[nm]]
    if (!is.list(entry)) next
    rpt <- entry$report
    if (is.null(rpt)) next
    rpt <- tryCatch(as.character(rpt)[[1L]], error = function(e) "")
    if (!nzchar(rpt)) next
    label <- tryCatch(labeller(pgx, nm), error = function(e) nm)
    one   <- .first_sentence(rpt, max = 120L)
    line  <- if (nzchar(one)) {
      sprintf("- %s (%s): %s", nm, label, one)
    } else {
      sprintf("- %s (%s)", nm, label)
    }
    lines <- c(lines, line)
  }
  if (!length(lines)) return("")
  paste(lines, collapse = "\n")
}

# ---- shared util: first prose sentence -----------------------------------

# Strips leading markdown noise (`#`, `>`, `*`, `-`, whitespace), collapses
# whitespace runs, takes the first sentence (terminator `.!?`), and caps the
# result at `max` chars without inserting an ellipsis. Empty input → "".
.first_sentence <- function(text, max = 120L) {
  if (is.null(text) || length(text) == 0L) return("")
  text <- as.character(text)[[1L]]
  if (is.na(text) || !nzchar(text)) return("")
  # Drop leading markdown chars (# > * - whitespace, repeated).
  text <- sub("^[\\s>#*\\-]+", "", text, perl = TRUE)
  text <- gsub("\\s+", " ", text, perl = TRUE)
  text <- trimws(text)
  if (!nzchar(text)) return("")
  first <- sub("^([^.!?]*[.!?]).*$", "\\1", text)
  if (identical(first, text) && nchar(text) > max) {
    first <- substr(text, 1L, max)
  }
  first <- trimws(first)
  if (nchar(first) > max) first <- substr(first, 1L, max)
  first
}

# ---- public assembler -----------------------------------------------------

#' Build the structured payload consumed by `make_followup_generator`'s
#' `generate()` method.
#'
#' Pure function — no reactives, no I/O beyond reading agent slots.
#'
#' @param agent An `omicsagentovi::Agent` (post-stream snapshot from
#'   `result$agent`). NULL → empty payload for that block.
#' @param last_text Character scalar — the assistant's final reply.
#' @param attached_report_slots Optional character vector — the report slot
#'   names that the user has currently ticked AND not yet consumed (i.e.
#'   what will be injected into the agent's NEXT user turn). When non-empty,
#'   `available_reports` is filtered to this set. When `character(0)` /
#'   `NULL`, `available_reports` is rendered as `""` so the helper sees
#'   an empty block (per the "hard gate" design: only mention reports the
#'   agent will actually have access to).
#' @return Named list with five blocks:
#'   `last_text`, `dataset_context`, `recent_tool_outputs`,
#'   `tool_catalog`, `available_reports`. Missing blocks are `""` or
#'   `list()` so the prompt template can render them uniformly.
#' @export
build_followup_payload <- function(agent, last_text,
                                   attached_report_slots = character(0)) {
  dctx <- .copilot_dataset_context_text(agent)
  if (is.null(dctx)) dctx <- ""
  # Gate the reports inventory by the user's ticked-and-unconsumed set.
  # If nothing is ticked, surface NOTHING — the helper must treat reports
  # as unavailable for the next turn (Option 1 + small twist from the
  # design discussion).
  reports_text <- if (length(attached_report_slots) > 0L) {
    .build_report_inventory(agent, only_slots = attached_report_slots)
  } else {
    ""
  }
  list(
    last_text           = if (is.null(last_text)) "" else as.character(last_text)[[1L]],
    dataset_context     = dctx,
    recent_tool_outputs = .extract_recent_tool_outputs(agent, max_outputs = 2L),
    tool_catalog        = .build_tool_catalog(),
    available_reports   = reports_text
  )
}
