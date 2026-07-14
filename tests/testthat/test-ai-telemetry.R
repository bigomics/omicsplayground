## Tests for the AI telemetry core module (components/modules/AiTelemetry.R).
## Covers pricing resolution + math, DB schema/idempotency, the turn parser,
## record() push sink, collect_chat() idempotent pull, and record_reports().

# AiTelemetry.R (sourced below into globalenv) uses %||%, whose canonical
# definition lives in components/app/R/utils/utils.R. source() defaults to the
# global env, so pull in that real definition rather than re-defining the
# operator here (project rule: never shim %||%).
.opg_root <- rprojroot::find_root(rprojroot::has_file("DESCRIPTION"))
source(file.path(.opg_root, "components", "app", "R", "utils", "utils.R"))
source(file.path(.opg_root, "components", "modules", "AiTelemetry.R"))

# A pricing YAML written into `dir` so tests don't depend on etc/.
.write_pricing_yaml <- function(dir) {
  path <- file.path(dir, "model_pricing.yml")
  writeLines(c(
    'version: "test-2026-06-23"',
    "models:",
    '  "*gpt-5.4-nano*":',
    "    input: 0.05",
    "    cached: 0.005",
    "    output: 0.40",
    '  "*":',
    "    input: 1.0",
    "    cached: 0.1",
    "    output: 2.0"
  ), path)
  path
}

# Build a faithful assistant payload via the real codec round-trip, so the
# parser is exercised against exactly what sessions.sqlite stores. ellmer
# normalises tokens to c(fresh, output, cached); reasoning is carried only in
# json$usage, so this also exercises the json fallback for reasoning_tokens.
.synthetic_payload <- function() {
  turn <- ellmer::Turn(
    role = "assistant",
    contents = list(
      ellmer::ContentText(text = "hello"),
      ellmer::ContentToolRequest(id = "t1", name = "query_de", arguments = list())
    ),
    tokens = c(13966, 37, 2000)
  )
  S7::prop(turn, "json") <- list(
    model = "gpt-5.4-nano-2026-03-17",
    usage = list(
      input_tokens_details  = list(cached_tokens = 2000),
      output_tokens_details = list(reasoning_tokens = 10)
    )
  )
  omicsagentovi::turn_to_payload(turn)
}

test_that("pricing resolves both the session and dated model to the same rates", {
  d <- tempfile("pdir"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  yml <- .write_pricing_yaml(d)
  a <- ai_telemetry_pricing("openai:gpt-5.4-nano", path = yml)
  b <- ai_telemetry_pricing("gpt-5.4-nano-2026-03-17", path = yml)
  expect_equal(a$input, 0.05)
  expect_equal(a[c("input", "cached", "output")], b[c("input", "cached", "output")])
  expect_equal(a$version, "test-2026-06-23")
})

test_that("catch-all glob is used only when nothing more specific matches", {
  d <- tempfile("pdir"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  yml <- .write_pricing_yaml(d)
  expect_equal(ai_telemetry_pricing("some-unknown-model", path = yml)$input, 1.0)
})

test_that("price matches a hand-computed reference", {
  pricing <- list(input = 0.05, cached = 0.005, output = 0.40, version = "x")
  usage <- list(input_fresh = 13966, input_cached = 2000, output = 37, reasoning = 10)
  ref <- 13966/1e6*0.05 + 2000/1e6*0.005 + (37 + 10)/1e6*0.40
  expect_equal(ai_telemetry_price(usage, "m", pricing), round(ref, 6))
})

test_that("price is NA (not 0) when no token detail is known", {
  pricing <- list(input = 0.05, cached = 0.005, output = 0.40, version = "x")
  usage <- list(input_fresh = NA, input_cached = NA, output = NA, reasoning = NA)
  expect_true(is.na(ai_telemetry_price(usage, "m", pricing)))
})

test_that("price is NA when pricing is unresolved", {
  expect_true(is.na(ai_telemetry_price(list(output = 5), "m", pricing = NULL)))
})

test_that("opening a fresh db creates the table and re-open is idempotent", {
  db <- tempfile("tel", fileext = ".sqlite"); on.exit(unlink(db))
  con <- ai_telemetry_open(db)
  expect_true("ai_usage_events" %in% DBI::dbListTables(con))
  DBI::dbDisconnect(con)
  con2 <- ai_telemetry_open(db)
  expect_true("ai_usage_events" %in% DBI::dbListTables(con2))
  DBI::dbDisconnect(con2)
})

test_that("parser extracts tokens, cache, reasoning and tool_count", {
  skip_if_not(requireNamespace("omicsai", quietly = TRUE))
  skip_if_not(requireNamespace("omicsagentovi", quietly = TRUE))
  p <- .ai_tel_parse_turn(.synthetic_payload(), session_model = "openai:gpt-5.4-nano")
  expect_equal(p$tool_count, 1L)
  expect_equal(p$model, "gpt-5.4-nano-2026-03-17")
  expect_equal(p$input_fresh, 13966L)
  expect_equal(p$output, 37L)
  expect_equal(p$input_cached, 2000L)  # recovered from json$usage
  expect_equal(p$reasoning, 10L)
})

test_that("parser never throws on malformed input", {
  p <- .ai_tel_parse_turn("{ not valid json", session_model = "m")
  expect_true(is.na(p$total))
  expect_equal(p$tool_count, 0L)
})

test_that("record() writes one event with the expected id and source", {
  d <- tempfile("recdir"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  db <- file.path(d, "ai_telemetry.sqlite")
  usage <- list(input_tokens_fresh = 100, input_tokens_cached = 0,
                output_tokens = 50, reasoning_tokens = 0)
  ai_telemetry_record("report", "sess-1", "u@x.com", "openai:gpt-5.4-nano",
                      usage, db_path = db)
  con <- DBI::dbConnect(RSQLite::SQLite(), db); on.exit(DBI::dbDisconnect(con), add = TRUE)
  rows <- DBI::dbGetQuery(con, "SELECT * FROM ai_usage_events")
  expect_equal(nrow(rows), 1L)
  expect_equal(rows$source, "report")
  expect_match(rows$event_id, "^report:sess-1:")  # source:session:timestamp
})

test_that("two one-shot events in the same session produce two rows (no collapse)", {
  d <- tempfile("recdir2"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  db <- file.path(d, "ai_telemetry.sqlite")
  usage <- list(input_tokens_fresh = 100, output_tokens = 50)
  t0 <- as.POSIXct(1750000000, origin = "1970-01-01", tz = "UTC")
  # Same source + session_id (a stable Shiny session token), distinct calls.
  ai_telemetry_record("card", "tok-1", "u@x.com", "openai:gpt-5.4-nano", usage,
                      created_at = t0,     db_path = db)
  ai_telemetry_record("card", "tok-1", "u@x.com", "openai:gpt-5.4-nano", usage,
                      created_at = t0 + 1, db_path = db)
  con <- DBI::dbConnect(RSQLite::SQLite(), db); on.exit(DBI::dbDisconnect(con), add = TRUE)
  n <- DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM ai_usage_events")$n
  expect_equal(n, 2L)  # would be 1 under the old session-only event_id
})

test_that("sessions db name + legacy migration helper", {
  expect_equal(.ai_tel_sessions_db_name("a@b.com"), "sessions-a@b.com.sqlite")
  expect_equal(.ai_tel_sessions_db_name(NA), "sessions.sqlite")

  d <- tempfile("chatdir"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  legacy <- file.path(d, "sessions.sqlite")
  writeLines("x", legacy)
  resolved <- .detect_legacy_sessions(d, "a@b.com")
  expect_equal(basename(resolved), "sessions-a@b.com.sqlite")
  expect_true(file.exists(resolved))
  expect_false(file.exists(legacy))
})

test_that("ai_telemetry_record_reports skips NULL-usage slots and records others", {
  d <- tempfile("repdir"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  db <- file.path(d, "ai_telemetry.sqlite")

  # usage mirrors what playbase stores: omicsai's native .omicsai_extract_usage()
  # field names + model. The sink reads these exact names, so the token columns
  # must populate — a renamed shape (total=/input=/output=) silently NA-fills.
  mk_usage <- function(total, fresh, cached, output) list(
    total_tokens = total, input_tokens = fresh + cached,
    input_tokens_fresh = fresh, input_tokens_cached = cached,
    output_tokens = output, reasoning_tokens = 0L, cache_hit_rate = 0,
    model = "openai:gpt-5.4-nano")

  pgx <- list(
    name = "test-dataset",
    ai = list(
      de = list(
        report     = "DE report",
        usage      = mk_usage(50L, 30L, 0L, 20L),
        created_at = 1750000000
      ),
      combined = list(
        report = "Combined report",
        usage  = NULL          # no usage — must be skipped
      ),
      pathways = list(
        report     = "Pathway report",
        usage      = mk_usage(80L, 55L, 5L, 20L),
        created_at = 1750000100
      )
    )
  )

  ids <- ai_telemetry_record_reports(pgx, user_email = NA_character_, db_path = db)
  con <- DBI::dbConnect(RSQLite::SQLite(), db); on.exit(DBI::dbDisconnect(con), add = TRUE)
  rows <- DBI::dbGetQuery(con, paste(
    "SELECT event_id, source, model, total, input_fresh, input_cached, output",
    "FROM ai_usage_events ORDER BY event_id"))

  # Only 2 slots have usage
  expect_equal(nrow(rows), 2L)
  expect_true(all(rows$source == "report"))
  # Both event_ids are dataset-keyed
  expect_true(all(grepl("^report:test-dataset:", rows$event_id)))
  # Regression guard: token columns must be populated, not NA (the field-name bug).
  de_row <- rows[grepl(":de:", rows$event_id), ]
  expect_equal(de_row$total, 50L)
  expect_equal(de_row$input_fresh, 30L)
  expect_equal(de_row$output, 20L)
  expect_equal(de_row$model, "openai:gpt-5.4-nano")
  expect_false(anyNA(rows$total))
})

test_that("ai_telemetry_record_reports is idempotent — second call inserts no new rows", {
  d <- tempfile("repdir2"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  db <- file.path(d, "ai_telemetry.sqlite")

  pgx <- list(
    name = "my-ds",
    ai = list(
      de = list(
        report     = "DE report",
        usage      = list(total = 50L, input = 30L, output = 20L,
                          model = "openai:gpt-5.4-nano"),
        created_at = 1750000000
      )
    )
  )

  ai_telemetry_record_reports(pgx, user_email = "u@x.com", db_path = db)
  ai_telemetry_record_reports(pgx, user_email = "u@x.com", db_path = db)

  con <- DBI::dbConnect(RSQLite::SQLite(), db); on.exit(DBI::dbDisconnect(con), add = TRUE)
  n <- DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM ai_usage_events")$n
  expect_equal(n, 1L)
})

test_that("ai_telemetry_record_reports tolerates NULL/empty pgx gracefully", {
  d <- tempfile("repdir3"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  db <- file.path(d, "ai_telemetry.sqlite")
  expect_silent(ai_telemetry_record_reports(NULL, NA_character_, db_path = db))
  expect_silent(ai_telemetry_record_reports(list(), NA_character_, db_path = db))
  expect_silent(ai_telemetry_record_reports(list(ai = list()), NA_character_, db_path = db))
})

test_that("collect_chat is idempotent over a synthetic chat store", {
  skip_if_not(requireNamespace("omicsai", quietly = TRUE))
  skip_if_not(requireNamespace("omicsagentovi", quietly = TRUE))
  d <- tempfile("chats"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  chat_db <- file.path(d, .ai_tel_sessions_db_name("u@x.com"))
  cc <- DBI::dbConnect(RSQLite::SQLite(), chat_db)
  DBI::dbExecute(cc, "CREATE TABLE sessions (session_id TEXT PRIMARY KEY, model TEXT)")
  DBI::dbExecute(cc, "CREATE TABLE turns (session_id TEXT, turn_idx INTEGER, role TEXT, payload_json TEXT, PRIMARY KEY (session_id, turn_idx))")
  DBI::dbExecute(cc, "INSERT INTO sessions VALUES ('s1', 'openai:gpt-5.4-nano')")
  pj <- jsonlite::toJSON(.synthetic_payload(), auto_unbox = TRUE, null = "null")
  DBI::dbExecute(cc, "INSERT INTO turns VALUES ('s1', 1, 'assistant', ?)", params = list(pj))
  DBI::dbExecute(cc, "INSERT INTO turns VALUES ('s1', 2, 'user', '{}')")
  DBI::dbDisconnect(cc)

  db <- file.path(d, "ai_telemetry.sqlite")
  n1 <- ai_telemetry_collect_chat(d, "u@x.com", db_path = db)
  n2 <- ai_telemetry_collect_chat(d, "u@x.com", db_path = db)
  expect_equal(n1, 1L)   # only the assistant turn
  expect_equal(n2, 0L)   # idempotent re-run
})
