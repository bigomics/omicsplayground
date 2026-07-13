## test-copilot-provider-args.R
##
## Unit tests for the copilot AI-provider / credential resolution helpers
## (.copilot_ai_provider, .copilot_agent_build_args) that thread the user's
## chosen provider + BYOK key into every Agent build / restore.
##
## These source the run controller from components/app_copilot/R and stub the
## app-global option accessors (getUserOption / get_ai_credentials / %||%) so
## the file loads standalone.

.board_dir <- if (dir.exists("components/app_copilot/R")) {
  "components/app_copilot/R"
} else {
  "../../components/app_copilot/R"
}

# ---- Standalone bootstrap: app-global helpers the module reads at runtime ----
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

# Minimal session double: userData is a shared env, exactly like Shiny's.
make_fake_session <- function(opts = list()) {
  ud <- new.env(parent = emptyenv())
  for (k in names(opts)) ud[[k]] <- opts[[k]]
  list(userData = ud)
}
getUserOption <- function(session, var, value) session$userData[[var]]
get_ai_credentials <- function(session) getUserOption(session, "ai_credentials")

source(file.path(.board_dir, "copilot_run_controller.R"), local = TRUE)

# ---------------------------------------------------------------------------
# .copilot_ai_provider
# ---------------------------------------------------------------------------
test_that(".copilot_ai_provider defaults to bigomics + NULL credentials", {
  s <- make_fake_session()
  sel <- .copilot_ai_provider(s)
  expect_equal(sel$provider, "bigomics")
  expect_null(sel$credentials)
})

test_that(".copilot_ai_provider surfaces a BYOK provider + credential closure", {
  cred <- local({ k <- "sk-123"; function() k })
  s <- make_fake_session(list(ai_provider = "anthropic", ai_credentials = cred))
  sel <- .copilot_ai_provider(s)
  expect_equal(sel$provider, "anthropic")
  expect_true(is.function(sel$credentials))
  expect_equal(sel$credentials(), "sk-123")
})

# ---------------------------------------------------------------------------
# .copilot_agent_build_args — BigOmics path (regression-critical)
# ---------------------------------------------------------------------------
test_that("BigOmics build args are byte-for-byte the tier path (no creds/model)", {
  s <- make_fake_session()  # provider defaults to bigomics
  for (tier in c("copilot-default", "copilot-deep")) {
    args <- .copilot_agent_build_args(s, tier)
    expect_identical(args, list(tier = tier))
    expect_false("credentials" %in% names(args))
    expect_false("model" %in% names(args))
    expect_false("provider" %in% names(args))
  }
})

# ---------------------------------------------------------------------------
# .copilot_agent_build_args — BYOK path
# ---------------------------------------------------------------------------
test_that("BYOK picks the per-tier menu model, provider-prefixed, with key", {
  cred <- local({ k <- "sk-anthropic"; function() k })
  s <- make_fake_session(list(
    ai_provider          = "anthropic",
    ai_credentials       = cred,
    llm_copilot_deep     = "claude-sonnet-4-6",
    llm_copilot_balanced = "claude-haiku-4-5"
  ))

  deep <- .copilot_agent_build_args(s, "copilot-deep")
  expect_equal(deep$model, "anthropic:claude-sonnet-4-6")
  expect_identical(deep$credentials, cred)
  expect_false("tier" %in% names(deep))

  bal <- .copilot_agent_build_args(s, "copilot-default")
  expect_equal(bal$model, "anthropic:claude-haiku-4-5")
  expect_identical(bal$credentials, cred)
})

test_that("BYOK with an empty menu falls back to the provider-aware tier path", {
  cred <- local({ k <- "sk-google"; function() k })
  s <- make_fake_session(list(
    ai_provider          = "google",
    ai_credentials       = cred,
    llm_copilot_deep     = "",   # empty -> fall back
    llm_copilot_balanced = NULL
  ))

  deep <- .copilot_agent_build_args(s, "copilot-deep")
  expect_equal(deep$tier, "copilot-deep")
  expect_equal(deep$provider, "google")
  expect_identical(deep$credentials, cred)
  expect_false("model" %in% names(deep))

  bal <- .copilot_agent_build_args(s, "copilot-default")
  expect_equal(bal$tier, "copilot-default")
  expect_equal(bal$provider, "google")
})
