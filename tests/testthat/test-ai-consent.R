## test-ai-consent.R
##
## Unit tests for the per-user AI data-sharing consent store
## (components/modules/AiConsent.R). The consent is the one AI setting that
## outlives a session, so these focus on the fail-closed reads: a missing,
## corrupt, or stale-version record must never read as an opt-in.

.modules_dir <- if (dir.exists("components/modules")) {
  "components/modules"
} else {
  "../../components/modules"
}

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

## getUserOption double, matching Shiny's session$userData semantics.
getUserOption <- function(session, var, value) session$userData[[var]]
make_fake_session <- function(opts = list()) {
  ud <- new.env(parent = emptyenv())
  for (k in names(opts)) ud[[k]] <- opts[[k]]
  list(userData = ud)
}

source(file.path(.modules_dir, "AiConsent.R"), local = TRUE)

# ---------------------------------------------------------------------------
# Round trip
# ---------------------------------------------------------------------------
test_that("consent round-trips through the user dir", {
  dir <- withr::local_tempdir()

  expect_false(load_ai_consent(dir))          # nothing written yet

  save_ai_consent(dir, TRUE)
  expect_true(file.exists(file.path(dir, AI_CONSENT_FILE)))
  expect_true(load_ai_consent(dir))

  save_ai_consent(dir, FALSE)
  expect_false(load_ai_consent(dir))
})

test_that("save records the version and a timestamp alongside the decision", {
  dir <- withr::local_tempdir()
  save_ai_consent(dir, TRUE)
  rec <- jsonlite::read_json(file.path(dir, AI_CONSENT_FILE), simplifyVector = TRUE)
  expect_true(rec$share_data)
  expect_equal(rec$version, AI_CONSENT_VERSION)
  expect_true(is.numeric(rec$consented_at) && rec$consented_at > 0)
})

# ---------------------------------------------------------------------------
# Fail-closed reads
# ---------------------------------------------------------------------------
test_that("load fails closed on missing, NULL and non-existent dirs", {
  expect_false(load_ai_consent(NULL))
  expect_false(load_ai_consent(""))
  expect_false(load_ai_consent(file.path(tempdir(), "no-such-user-dir")))
})

test_that("a corrupt record reads as no consent", {
  dir <- withr::local_tempdir()
  writeLines("{not json at all", file.path(dir, AI_CONSENT_FILE))
  expect_false(load_ai_consent(dir))
})

test_that("a record from an older consent version reads as no consent", {
  ## The wording changed materially, so a previously recorded opt-in is no
  ## longer informed — the user must be asked again.
  dir <- withr::local_tempdir()
  jsonlite::write_json(
    list(share_data = TRUE, version = AI_CONSENT_VERSION - 1L, consented_at = 1),
    file.path(dir, AI_CONSENT_FILE), auto_unbox = TRUE
  )
  expect_false(load_ai_consent(dir))
})

test_that("a record with no version field reads as no consent", {
  dir <- withr::local_tempdir()
  jsonlite::write_json(list(share_data = TRUE), file.path(dir, AI_CONSENT_FILE),
                       auto_unbox = TRUE)
  expect_false(load_ai_consent(dir))
})

# ---------------------------------------------------------------------------
# Writes to an unwritable / absent dir are a no-op, not an error
# ---------------------------------------------------------------------------
test_that("saving to a non-existent user dir is a silent no-op", {
  ## No-auth deployments have no user_dir; the consent simply stays
  ## session-scoped rather than blowing up the settings observer.
  expect_silent(save_ai_consent(NULL, TRUE))
  expect_silent(save_ai_consent(file.path(tempdir(), "nope"), TRUE))
})

# ---------------------------------------------------------------------------
# Session accessor
# ---------------------------------------------------------------------------
test_that("get_ai_data_consent reads the session copy and defaults to FALSE", {
  expect_false(get_ai_data_consent(make_fake_session()))
  expect_false(get_ai_data_consent(make_fake_session(list(ai_share_data = NULL))))
  expect_false(get_ai_data_consent(make_fake_session(list(ai_share_data = NA))))
  expect_true(get_ai_data_consent(make_fake_session(list(ai_share_data = TRUE))))
})
