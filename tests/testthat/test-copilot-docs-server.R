## test-copilot-docs-server.R

.board_dir <- if (dir.exists("components/app_copilot/R")) {
  "components/app_copilot/R"
} else {
  "../../components/app_copilot/R"
}

source(file.path(.board_dir, "CopilotDocsServer.R"), local = TRUE)

# Tests for the Copilot document upload gate.
#
# Regression origin: a customer uploaded a .json and a 6 MB .txt. Both were
# listed in the document panel and neither ever reached the agent, so the
# assistant answered as if only the PDF existed and the user concluded the
# chatbot could not read .txt at all. The invariant these tests protect is
# that the panel lists exactly what the agent can ingest.

test_that(".copilot_filter_uploads accepts the three supported extensions", {
  v <- .copilot_filter_uploads(
    c("paper.pdf", "notes.txt", "readme.md"),
    c(1000, 1000, 1000)
  )
  expect_true(all(v$ok))
  expect_true(all(is.na(v$reason)))
})

test_that(".copilot_filter_uploads is case-insensitive on the extension", {
  v <- .copilot_filter_uploads(c("PAPER.PDF", "Notes.TxT"), c(10, 10))
  expect_true(all(v$ok))
})

test_that(".copilot_filter_uploads rejects the .json that started this", {
  v <- .copilot_filter_uploads(
    c("corum_humanComplexes.json", "paper.pdf"),
    c(15600 * 1024, 3915 * 1024)
  )
  expect_equal(v$ok, c(FALSE, TRUE))
  expect_match(v$reason[[1L]], "\\.json is not supported")
  expect_match(v$reason[[1L]], "\\.pdf")
  expect_true(is.na(v$reason[[2L]]))
})

test_that(".copilot_filter_uploads names the extension it saw", {
  v <- .copilot_filter_uploads(c("data.csv", "sheet.xlsx"), c(10, 10))
  expect_false(any(v$ok))
  expect_match(v$reason[[1L]], "\\.csv")
  expect_match(v$reason[[2L]], "\\.xlsx")
})

test_that(".copilot_filter_uploads handles a file with no extension", {
  v <- .copilot_filter_uploads("README", 100)
  expect_false(v$ok)
  expect_match(v$reason, "no extension")
})

test_that(".copilot_filter_uploads enforces the size ceiling", {
  v <- .copilot_filter_uploads(
    c("small.txt", "huge.txt"),
    c(.COPILOT_DOC_MAX_BYTES - 1, .COPILOT_DOC_MAX_BYTES + 1)
  )
  expect_equal(v$ok, c(TRUE, FALSE))
  expect_match(v$reason[[2L]], "too large")
  expect_match(v$reason[[2L]], "20 MB")
})

test_that("the 6 MB table from the customer report is now accepted", {
  # It was rejected by the agent's old 500 KB cap; the cap is 20 MB now and
  # the UI must not re-impose the old limit.
  v <- .copilot_filter_uploads("corum_humanComplexes.txt", 6122 * 1024)
  expect_true(v$ok)
})

test_that(".copilot_filter_uploads rejects empty files", {
  v <- .copilot_filter_uploads("empty.txt", 0)
  expect_false(v$ok)
  expect_match(v$reason, "empty")
})

test_that(".copilot_filter_uploads reports extension before size", {
  # A huge .json should be told it is the wrong type, not the wrong size —
  # shrinking it would not help.
  v <- .copilot_filter_uploads("big.json", .COPILOT_DOC_MAX_BYTES * 2)
  expect_false(v$ok)
  expect_match(v$reason, "not supported")
})

test_that(".copilot_filter_uploads tolerates a missing size", {
  v <- .copilot_filter_uploads("notes.txt", NA)
  expect_true(v$ok)
})

test_that(".copilot_filter_uploads handles an empty upload set", {
  v <- .copilot_filter_uploads(character(0), numeric(0))
  expect_length(v$ok, 0L)
  expect_length(v$reason, 0L)
})

test_that("the listing pattern matches exactly the agent's ingest set", {
  # .COPILOT_DOC_PATTERN is the panel's view of the directory and must agree
  # with omicsagentovi's bootstrap glob, or the panel shows files the agent
  # cannot see (the original bug) or hides files it can.
  supported   <- c("a.pdf", "b.txt", "c.md", "D.PDF", "E.Md")
  unsupported <- c("x.json", "y.csv", "z.xlsx", "no_ext", "w.pdf.bak")

  expect_true(all(grepl(.COPILOT_DOC_PATTERN, supported, ignore.case = TRUE)))
  expect_false(any(grepl(.COPILOT_DOC_PATTERN, unsupported, ignore.case = TRUE)))

  # And the pattern agrees with the extension list the gate uses.
  expect_setequal(.COPILOT_DOC_EXTS, c("pdf", "txt", "md"))
})

test_that("the UI cap agrees with the agent's own limit", {
  # The whole defect class here is the UI and the agent disagreeing about
  # what is ingestible. DESCRIPTION pins omicsagentovi (>= 0.5.90) but
  # Suggests is not enforced at runtime, so assert the two constants
  # directly whenever a new enough agent is actually installed.
  testthat::skip_if_not_installed("omicsagentovi", "0.5.90")
  backend <- tryCatch(omicsagentovi:::.OVI_MAX_DOC_CHARS,
                      error = function(e) NULL)
  testthat::skip_if(is.null(backend), "agent does not expose its doc limit")

  # Bytes here approximate the agent's extracted-character limit; the UI
  # must never allow more than the agent will accept.
  expect_lte(.COPILOT_DOC_MAX_BYTES, backend)
})

test_that("the docs listing hides an unsupported file already on disk", {
  # Defends the panel against files staged by an older build or dropped in
  # by hand: what is listed must be what the agent globs.
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  writeLines("body", file.path(dir, "paper.txt"))
  writeLines("{}",   file.path(dir, "corum.json"))

  listed <- list.files(dir, pattern = .COPILOT_DOC_PATTERN, ignore.case = TRUE)
  expect_equal(listed, "paper.txt")
})
