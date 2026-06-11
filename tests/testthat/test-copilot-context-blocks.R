## test-copilot-context-blocks.R

.board_dir <- if (dir.exists("components/app_copilot/R")) {
  "components/app_copilot/R"
} else {
  "../../components/app_copilot/R"
}
.module_dir <- if (dir.exists("components/modules")) {
  "components/modules"
} else {
  "../../components/modules"
}

source(file.path(.board_dir, "copilot_logger.R"), local = TRUE)
source(file.path(.module_dir, "AiReports.R"), local = TRUE)
source(file.path(.board_dir, "CopilotReportsServer.R"), local = TRUE)
source(file.path(.board_dir, "copilot_context_blocks.R"), local = TRUE)

.pgx_with_ai <- function(ai) {
  structure(
    list(
      name = "demo",
      date = "2026-06-11",
      X = matrix(1, nrow = 2, ncol = 2),
      ai = ai
    ),
    class = "pgx"
  )
}

test_that("ai report context uses combined first when available", {
  pgx <- .pgx_with_ai(list(
    de = list(report = "DE report"),
    combined = list(report = "Combined report")
  ))

  out <- .copilot_ai_report_context(pgx)

  expect_match(out, "Summary \\[combined\\]")
  expect_match(out, "Combined report", fixed = TRUE)
  expect_false(grepl("DE report", out, fixed = TRUE))
})

test_that("ai report context falls back to module snippets", {
  pgx <- .pgx_with_ai(list(
    pathways = list(report = "Pathway report"),
    de = list(report = "DE report")
  ))

  out <- .copilot_ai_report_context(pgx)

  expect_match(out, "Differential Expression \\[de\\]")
  expect_match(out, "Enrichment \\[pathways\\]")
  expect_match(out, "DE report", fixed = TRUE)
  expect_match(out, "Pathway report", fixed = TRUE)
})

test_that("ai report context returns NULL without usable pgx ai reports", {
  expect_null(.copilot_ai_report_context(list()))
  expect_null(.copilot_ai_report_context(.pgx_with_ai(list(
    de = list(report = "")
  ))))
})

test_that("ai report context enforces max_chars budget", {
  pgx <- .pgx_with_ai(list(
    combined = list(report = paste(rep("x", 500), collapse = ""))
  ))

  out <- .copilot_ai_report_context(pgx, max_chars = 160L)

  expect_lte(nchar(out, type = "chars"), 160L)
  expect_match(out, "truncated", fixed = TRUE)
})

test_that("explicit report slot selection controls included reports", {
  pgx <- .pgx_with_ai(list(
    combined = list(report = "Combined report"),
    de = list(report = "DE report")
  ))

  out <- .copilot_ai_report_context(pgx, slots = "de")

  expect_match(out, "DE report", fixed = TRUE)
  expect_false(grepl("Combined report", out, fixed = TRUE))
})
