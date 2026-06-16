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

test_that("user docs context concatenates selected files as labeled sections", {
  tmp <- tempfile("copilot_docs_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  writeLines("alpha body", file.path(tmp, "alpha.txt"))
  writeLines(c("# Beta", "beta body"), file.path(tmp, "beta.md"))

  built <- .copilot_user_docs_context(tmp, c("alpha.txt", "beta.md"))

  expect_match(built$text, "## alpha.txt", fixed = TRUE)
  expect_match(built$text, "alpha body", fixed = TRUE)
  expect_match(built$text, "## beta.md", fixed = TRUE)
  expect_match(built$text, "beta body", fixed = TRUE)
  expect_setequal(built$docs, c("alpha.txt", "beta.md"))
})

test_that("user docs context applies no character cap", {
  tmp <- tempfile("copilot_docs_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  big <- paste(rep("x", 50000L), collapse = "")
  writeLines(big, file.path(tmp, "big.txt"))

  built <- .copilot_user_docs_context(tmp, "big.txt")

  expect_gt(nchar(built$text, type = "chars"), 50000L)
  expect_false(grepl("[truncated]", built$text, fixed = TRUE))
})

test_that("user docs context skips missing or empty files", {
  tmp <- tempfile("copilot_docs_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  writeLines("only one", file.path(tmp, "ok.txt"))

  built <- .copilot_user_docs_context(tmp, c("ok.txt", "missing.txt"))

  expect_match(built$text, "ok.txt", fixed = TRUE)
  expect_false(grepl("missing.txt", built$text, fixed = TRUE))
  expect_setequal(built$docs, "ok.txt")
})

test_that("user docs context returns NULL text when nothing usable", {
  tmp <- tempfile("copilot_docs_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  built <- .copilot_user_docs_context(tmp, c("nope.txt"))

  expect_null(built$text)
  expect_length(built$docs, 0L)
})

# ===========================================================================
# .copilot_dataset_context_text — extracted helper shared by the
# current_dataset provider and the follow-up payload builder.
# ===========================================================================

.make_dataset_agent <- function(name, organism, samples, genes, contrasts,
                                description) {
  CtxC <- S7::new_class("DatasetCtxStub",
                        properties = list(
                          pgx          = S7::class_any,
                          dataset_name = S7::class_character
                        ))
  AgC  <- S7::new_class("DatasetAgentStub",
                        properties = list(context = CtxC))
  pgx  <- structure(
    list(
      samples     = matrix(0, nrow = samples, ncol = 1L),
      X           = matrix(0, nrow = genes,   ncol = 1L),
      contrasts   = setNames(vector("list", length(contrasts)), contrasts),
      organism    = organism,
      description = description
    ),
    class = "pgx"
  )
  AgC(context = CtxC(pgx = pgx, dataset_name = name))
}

test_that(".copilot_dataset_context_text builds labeled lines for a full pgx", {
  agent <- .make_dataset_agent(
    name        = "example-data",
    organism    = "Homo sapiens",
    samples     = 12L,
    genes       = 18000L,
    contrasts   = c("act48h_vs_notact", "ctrl_vs_treated"),
    description = "Activation timecourse."
  )
  out <- .copilot_dataset_context_text(agent)
  expect_match(out, "name: example-data",             fixed = TRUE)
  expect_match(out, "organism: Homo sapiens",         fixed = TRUE)
  expect_match(out, "samples: 12",                    fixed = TRUE)
  expect_match(out, "genes: 18000",                   fixed = TRUE)
  expect_match(out, "act48h_vs_notact, ctrl_vs_treated", fixed = TRUE)
  expect_match(out, "description: Activation",        fixed = TRUE)
})

test_that(".copilot_dataset_context_text returns NULL for NULL agent / no pgx", {
  expect_null(.copilot_dataset_context_text(NULL))
})
