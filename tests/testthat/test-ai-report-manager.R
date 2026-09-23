##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## Tests for the app-scope AI report manager (components/modules/AiReportManager.R).
## Covers run-key identity, the subscriber registry, sidecar addressing and the
## merge semantics of ai_report_commit_to_disk().

.opg_root <- rprojroot::find_root(rprojroot::has_file("DESCRIPTION"))

## The manager logs through info(), which lives in components/utils/utils.R -
## a file that also configures waiter at source time. Nothing here asserts on
## log output, so a quiet stand-in keeps the test process free of that.
info <- function(..., type = "INFO") invisible(NULL)

source(file.path(.opg_root, "components", "modules", "AiReportManager.R"),
       local = TRUE)

## A run registered by hand: the real ai_report_run_start() dispatches to mirai
## workers, and none of the registry behaviour under test needs that.
.fake_run <- function(key) {
  run <- new.env(parent = emptyenv())
  run$running     <- TRUE
  run$done        <- 2L
  run$total       <- 5L
  run$failed      <- 1L
  run$slot        <- "de"
  run$started     <- Sys.time()
  run$subscribers <- list()
  run$next_sub    <- 0L
  .ai_report_runs[[key]] <- run
  withr::defer(
    suppressWarnings(rm(list = key, envir = .ai_report_runs)),
    envir = parent.frame()
  )
  run
}

.write_pgx <- function(dir, ai) {
  path <- file.path(dir, "testset.pgx")
  pgx <- list(name = "testset", X = matrix(0, 2, 2))
  if (!missing(ai)) pgx$ai <- ai
  suppressMessages(capture.output(playbase::pgx.save(pgx, file = path)))
  path
}

test_that("run key separates users and models on the same dataset", {
  path <- file.path(tempdir(), "shared.pgx")
  base <- ai_report_run_key(path, "set|10x4|2026-01-01", "ann@x.io", "gpt-5")

  expect_identical(
    base,
    ai_report_run_key(path, "set|10x4|2026-01-01", "ann@x.io", "gpt-5")
  )
  expect_false(base ==
    ai_report_run_key(path, "set|10x4|2026-01-01", "bob@x.io", "gpt-5"))
  expect_false(base ==
    ai_report_run_key(path, "set|10x4|2026-01-01", "ann@x.io", "claude-4"))
  expect_false(base ==
    ai_report_run_key(path, "set|11x4|2026-01-02", "ann@x.io", "gpt-5"))
  expect_false(base ==
    ai_report_run_key(file.path(tempdir(), "other.pgx"),
                      "set|10x4|2026-01-01", "ann@x.io", "gpt-5"))
})

test_that("run key tolerates a missing user or model", {
  path <- file.path(tempdir(), "shared.pgx")
  anon <- ai_report_run_key(path, "tok", NULL, "gpt-5")

  expect_type(anon, "character")
  expect_length(anon, 1L)
  expect_false(anon == ai_report_run_key(path, "tok", "ann@x.io", "gpt-5"))
})

test_that("each subscription gets its own id and its own callbacks", {
  run <- .fake_run("key-subs")
  seen <- character(0)

  first <- ai_report_run_subscribe("key-subs",
    on_done = function(result) seen <<- c(seen, "first"))
  second <- ai_report_run_subscribe("key-subs",
    on_done = function(result) seen <<- c(seen, "second"))

  expect_false(first == second)
  expect_length(run$subscribers, 2L)

  ## Two tabs of one Shiny session used to share a subscriber id, so the
  ## second subscription silently replaced the first one's callbacks.
  .ai_report_notify(run, "on_done", list(ai = NULL))
  expect_setequal(seen, c("first", "second"))
})

test_that("a dead subscriber does not stop the others being notified", {
  run <- .fake_run("key-throw")
  seen <- character(0)

  ai_report_run_subscribe("key-throw",
    on_progress = function(...) stop("session gone"))
  ai_report_run_subscribe("key-throw",
    on_progress = function(...) seen <<- c(seen, "alive"))

  expect_silent(.ai_report_notify(run, "on_progress", 1L, 2L, "de", TRUE))
  expect_identical(seen, "alive")
})

test_that("subscribe returns NULL when there is no run to join", {
  expect_null(ai_report_run_subscribe("key-absent"))
  expect_false(ai_report_run_unsubscribe("key-absent", "sub1"))
})

test_that("unsubscribe drops only the given subscriber", {
  run <- .fake_run("key-unsub")
  seen <- character(0)

  first <- ai_report_run_subscribe("key-unsub",
    on_done = function(result) seen <<- c(seen, "first"))
  second <- ai_report_run_subscribe("key-unsub",
    on_done = function(result) seen <<- c(seen, "second"))

  expect_true(ai_report_run_unsubscribe("key-unsub", first))
  expect_length(run$subscribers, 1L)

  .ai_report_notify(run, "on_done", list(ai = NULL))
  expect_identical(seen, "second")

  ## An unknown id is harmless: a session may end after the run unregistered.
  expect_true(ai_report_run_unsubscribe("key-unsub", "sub99"))
  expect_length(run$subscribers, 1L)
})

test_that("run status reports progress and disappears with the run", {
  run <- .fake_run("key-status")

  status <- ai_report_run_status("key-status")
  expect_identical(status$running, TRUE)
  expect_identical(status$done, 2L)
  expect_identical(status$total, 5L)
  expect_identical(status$failed, 1L)
  expect_identical(status$slot, "de")
  expect_s3_class(status$started, "POSIXct")

  expect_true(ai_report_run_active("key-status"))
  run$running <- FALSE
  expect_false(ai_report_run_active("key-status"))
  expect_null(ai_report_run_status("key-absent"))
})

test_that("sidecar directory is keyed on the dataset token, not the path alone", {
  path <- file.path(tempdir(), "sets", "testset.pgx")
  v1 <- ai_report_sidecar_dir(path, "testset|100x8|2026-01-01")
  v2 <- ai_report_sidecar_dir(path, "testset|120x8|2026-02-01")

  expect_identical(v1, ai_report_sidecar_dir(path, "testset|100x8|2026-01-01"))
  expect_false(v1 == v2)
  expect_identical(dirname(dirname(v1)), dirname(path))
  expect_identical(basename(dirname(v1)), ".ai_reports")
  expect_true(startsWith(basename(v1), "testset-"))
})

test_that("sidecars of two dataset versions do not overwrite each other", {
  dir <- withr::local_tempdir()
  path <- file.path(dir, "testset.pgx")

  ai_report_write_sidecar(path, "testset|100x8|2026-01-01", "de", "old report")
  ai_report_write_sidecar(path, "testset|120x8|2026-02-01", "de", "new report")

  old <- file.path(ai_report_sidecar_dir(path, "testset|100x8|2026-01-01"), "de.md")
  new <- file.path(ai_report_sidecar_dir(path, "testset|120x8|2026-02-01"), "de.md")
  expect_identical(readLines(old), "old report")
  expect_identical(readLines(new), "new report")

  meta <- jsonlite::fromJSON(sub("[.]md$", ".json", new))
  expect_identical(meta$slot, "de")
  expect_identical(meta$chars, nchar("new report"))
})

test_that("commit merges into the reports already on disk", {
  dir <- withr::local_tempdir()
  path <- .write_pgx(dir, ai = list(
    de = list(report = "old de"),
    clustering = list(report = "old clustering")
  ))

  expect_true(ai_report_commit_to_disk(path, list(
    de = list(report = "new de"),
    pathways = list(report = "new pathways")
  )))

  ai <- playbase::pgx.load(path)$ai
  expect_identical(ai$de$report, "new de")
  expect_identical(ai$pathways$report, "new pathways")
  expect_identical(ai$clustering$report, "old clustering")
})

test_that("commit never deletes a slot it was handed as NULL", {
  dir <- withr::local_tempdir()
  path <- .write_pgx(dir, ai = list(de = list(report = "old de")))

  ## `existing[[slot]] <- NULL` would drop the report that is already saved.
  expect_true(ai_report_commit_to_disk(path,
    list(de = NULL, pathways = list(report = "new pathways"))))

  ai <- playbase::pgx.load(path)$ai
  expect_identical(ai$de$report, "old de")
  expect_identical(ai$pathways$report, "new pathways")
})

test_that("commit reports FALSE when there is nothing to write", {
  dir <- withr::local_tempdir()
  path <- .write_pgx(dir, ai = list(de = list(report = "old de")))
  before <- file.mtime(path)

  expect_false(ai_report_commit_to_disk(path, NULL))
  expect_false(ai_report_commit_to_disk(path, list()))
  expect_false(ai_report_commit_to_disk(path, list(de = NULL)))
  expect_false(ai_report_commit_to_disk(path, "not a list"))
  expect_false(ai_report_commit_to_disk(file.path(dir, "absent.pgx"),
                                        list(de = list(report = "x"))))

  expect_identical(file.mtime(path), before)
  expect_identical(playbase::pgx.load(path)$ai$de$report, "old de")
})

test_that("commit creates the ai slot when the pgx has none or a broken one", {
  dir <- withr::local_tempdir()
  missing_ai <- .write_pgx(dir)
  broken_ai <- file.path(dir, "broken.pgx")
  file.copy(.write_pgx(withr::local_tempdir(), ai = "not a list"), broken_ai)

  expect_true(ai_report_commit_to_disk(missing_ai,
    list(de = list(report = "fresh"))))
  expect_identical(playbase::pgx.load(missing_ai)$ai$de$report, "fresh")

  expect_true(ai_report_commit_to_disk(broken_ai,
    list(de = list(report = "fresh"))))
  ai <- playbase::pgx.load(broken_ai)$ai
  expect_type(ai, "list")
  expect_identical(ai$de$report, "fresh")
})
