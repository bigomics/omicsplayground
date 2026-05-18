## test-copilot-bindings.R
##
## Tier 1 pure unit tests for build_run_bindings().
## Matches the test plan in .active_plans/refactor_copilot/bindings/specs.md §"Tier 1".
## No live LLM, no Shiny server loop required.
##
## Requires omicsagentovi 0.4.x (agent_refactor branch) loaded from source.
## The test suite skips cleanly if only 0.3.x is available.

# Load omicsagentovi from source (0.4.x S7 API) if the installed version is too old.
.ovi_source_dir <- Sys.getenv("OVI_SOURCE_DIR", "/Users/santiago/projects/BigOmics/omicsagentovi")
if (packageVersion("omicsagentovi") < "0.4.0" && dir.exists(.ovi_source_dir)) {
  pkgload::load_all(.ovi_source_dir, quiet = TRUE)
}

if (packageVersion("omicsagentovi") < "0.4.0") {
  skip("omicsagentovi >= 0.4.0 required (agent_refactor branch). Set OVI_SOURCE_DIR to the source path.")
}

# Source the factory under test.
# testthat sets the working directory to tests/testthat/, so we step up two
# levels to the repo root. devtools::test() also anchors at the repo root,
# so we check both.
.bindings_path <- if (file.exists("components/board.copilot/R/copilot_bindings.R")) {
  "components/board.copilot/R/copilot_bindings.R"
} else {
  "../../components/board.copilot/R/copilot_bindings.R"
}
source(.bindings_path, local = TRUE)

library(omicsagentovi)

# ── helpers ──────────────────────────────────────────────────────────────────

# Source the plot render helpers used by the new plot_callback implementation.
.render_path <- if (file.exists("components/board.copilot/R/copilot_plot_render.R")) {
  "components/board.copilot/R/copilot_plot_render.R"
} else {
  "../../components/board.copilot/R/copilot_plot_render.R"
}
source(.render_path, local = TRUE)

# A minimal evidence_api stub that records calls.
# post-Phase-5: $append_artifact(record) receives a single list arg.
make_evidence_stub <- function() {
  calls <- list()
  list(
    append_artifact = function(record) {
      calls[[length(calls) + 1L]] <<- record
      invisible(NULL)
    },
    calls = function() calls
  )
}

# A mock Shiny session with wrapFunction support.
make_mock_session <- function() {
  env <- new.env(parent = emptyenv())
  env$wrapFunction <- function(f) f   # identity — sufficient for unit tests
  env
}

# ── tests ─────────────────────────────────────────────────────────────────────

test_that("factory returns RunBindings S7 object", {
  result <- build_run_bindings(
    session          = make_mock_session(),
    evidence_api     = make_evidence_stub(),
    pgx_loaded_event = shiny::reactiveVal(NULL)
  )
  expect_true(S7::S7_inherits(result, omicsagentovi::RunBindings))
})

test_that("NULL session -> notification_sink is NULL", {
  result <- build_run_bindings(
    session          = NULL,
    evidence_api     = make_evidence_stub(),
    pgx_loaded_event = shiny::reactiveVal(NULL)
  )
  expect_null(result@notification_sink)
})

test_that("NULL session -> plot_callback is wrapped but functional (evidence_api present)", {
  # When session is NULL but evidence_api is provided, the impl closure is used directly.
  result <- build_run_bindings(
    session          = NULL,
    evidence_api     = make_evidence_stub(),
    pgx_loaded_event = NULL
  )
  # plot_callback should still be a function (no session$wrapFunction path taken)
  expect_true(is.function(result@plot_callback))
})

test_that("NULL evidence_api -> plot_callback is NULL", {
  result <- build_run_bindings(
    session          = make_mock_session(),
    evidence_api     = NULL,
    pgx_loaded_event = NULL
  )
  expect_null(result@plot_callback)
})

test_that("valid evidence_api -> plot_callback is a function", {
  result <- build_run_bindings(
    session          = make_mock_session(),
    evidence_api     = make_evidence_stub(),
    pgx_loaded_event = NULL
  )
  expect_true(is.function(result@plot_callback))
})

test_that("docs_dir NULL normalised to character(0)", {
  result <- build_run_bindings(
    session  = NULL,
    evidence_api = NULL,
    docs_dir = NULL
  )
  expect_identical(result@docs_dir, character(0))
})

test_that('docs_dir "" normalised to character(0)', {
  result <- build_run_bindings(
    session      = NULL,
    evidence_api = NULL,
    docs_dir     = ""
  )
  expect_identical(result@docs_dir, character(0))
})

test_that("docs_dir valid path passes through unchanged", {
  result <- build_run_bindings(
    session      = NULL,
    evidence_api = NULL,
    docs_dir     = "/tmp/docs"
  )
  expect_identical(result@docs_dir, "/tmp/docs")
})

test_that('notification_sink "pgx_loaded" writes payload to pgx_loaded_event reactiveVal', {
  # Use shiny::isolate to read reactiveVal outside a reactive context in tests.
  ev <- shiny::reactiveVal(NULL)
  result <- build_run_bindings(
    session          = make_mock_session(),
    evidence_api     = NULL,
    pgx_loaded_event = ev
  )
  payload <- list(
    pgx          = list(name = "mock"),
    dataset_name = "mock-ds",
    name_arg     = "/path/mock.pgx",
    data_dir     = "/path"
  )
  result@notification_sink("pgx_loaded", payload)
  received <- shiny::isolate(ev())
  expect_identical(received$dataset_name, "mock-ds")
  expect_identical(received$name_arg, "/path/mock.pgx")
})

test_that("notification_sink unknown level is a no-op (pgx_loaded_event unchanged)", {
  ev <- shiny::reactiveVal(NULL)
  result <- build_run_bindings(
    session          = make_mock_session(),
    evidence_api     = NULL,
    pgx_loaded_event = ev
  )
  result@notification_sink("some_unknown_level", list(foo = "bar"))
  expect_null(shiny::isolate(ev()))
})

test_that("pgx_loaded_event non-function non-reactiveVal -> error at factory call", {
  expect_error(
    build_run_bindings(
      session          = make_mock_session(),
      evidence_api     = NULL,
      pgx_loaded_event = "not_a_function_or_reactiveVal"
    ),
    regexp = "pgx_loaded_event must be"
  )
})

test_that("evidence_api missing append_artifact -> error at factory call", {
  expect_error(
    build_run_bindings(
      session      = make_mock_session(),
      evidence_api = list(something_else = function() NULL)
    ),
    regexp = "evidence_api must have"
  )
})

test_that("progress_callback is always a no-op function (never NULL)", {
  result <- build_run_bindings(session = NULL, evidence_api = NULL)
  expect_true(is.function(result@progress_callback))
  # Must not error and must return invisible(NULL)
  out <- result@progress_callback(list(kind = "tool_start", tool = "show_plot"))
  expect_null(out)
})

test_that("plot_callback propagates a well-formed record to evidence_api$append_artifact", {
  # Build a ggplot stub and override copilot_build_plot in the env where
  # copilot_bindings.R was sourced (same local env due to source(..., local = TRUE)).
  stub_plot <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
                 ggplot2::geom_point()
  old_build <- copilot_build_plot           # saved from the sourced env
  copilot_build_plot <<- function(...) stub_plot   # patch in same env
  on.exit(copilot_build_plot <<- old_build)        # restore

  stub   <- make_evidence_stub()
  result <- build_run_bindings(
    session          = NULL,
    evidence_api     = stub,
    pgx_loaded_event = NULL
  )

  # Invoke the callback directly; session = NULL means no wrapFunction path taken.
  result@plot_callback(
    pgx       = list(),
    plot_type = "pca",
    args      = list(contrast = "A_vs_B"),
    artifact  = NULL
  )

  recorded <- stub$calls()
  expect_equal(length(recorded), 1L)
  rec <- recorded[[1L]]
  expect_true(is.list(rec))
  expect_equal(rec$kind, "ggplot")
  expect_equal(rec$plot_type, "pca")
  expect_equal(rec$args$contrast, "A_vs_B")
  expect_null(rec$artifact)
  # Prerendered path must point to an existing PNG file
  expect_false(is.null(rec$prerendered_path))
  expect_true(file.exists(rec$prerendered_path))
  copilot_prerender_cleanup(rec$prerendered_path)
})
