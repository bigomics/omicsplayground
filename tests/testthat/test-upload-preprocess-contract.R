# Canonical bulk-upload preprocessing contract tests.
#
# These tests cover staged previews, pristine source reuse, and alignment.
# They exercise the installed playbase boundary used by the application.

.opg_test_root <- normalizePath(
  file.path(testthat::test_path(), "..", ".."),
  mustWork = TRUE
)
source(
  file.path(
    .opg_test_root,
    "components",
    "board.upload",
    "R",
    "upload_preprocess.R"
  ),
  local = TRUE
)

testthat::test_that("upload options use only the canonical vocabulary", {
  counts <- matrix(
    seq_len(24),
    nrow = 6,
    dimnames = list(
      c("gx:a", "gx:b", "gx:c", "px:a", "px:b", "px:c"),
      paste0("s", seq_len(4))
    )
  )
  options <- .opg_upload_preprocess_options(
    counts,
    datatype = "multi-omics",
    normalize = TRUE,
    norm_method = "multiomics",
    dedup = "off"
  )

  testthat::expect_setequal(
    names(options),
    c(
      "input_space",
      "output_space",
      "layers",
      "is_npx",
      "zero_as_na",
      "filter_missing",
      "filter_threshold",
      "impute",
      "impute_method",
      "impute_args",
      "normalize",
      "norm_method",
      "normalize_args",
      "remove_outliers",
      "outlier_threshold",
      "outlier_methods",
      "dedup",
      "batch_correct",
      "batch_method",
      "batch",
      "target",
      "batch_args",
      "max_features"
    )
  )
  testthat::expect_identical(
    options$norm_method,
    c(gx = "CPM", px = "maxMedian")
  )
  testthat::expect_identical(options$layers, c(rep("gx", 3), rep("px", 3)))
  testthat::expect_false(.opg_preview_options(options, "impute")$normalize)
  testthat::expect_false(
    .opg_preview_options(options, "normalize")$remove_outliers
  )
  testthat::expect_false(
    .opg_preview_options(options, "outliers")$batch_correct
  )
})

testthat::test_that("full upload preprocessing reruns without a result ratchet", {
  testthat::skip_if_not_installed("limma")
  testthat::expect_true(all(
    c(
      "pgx.preprocess",
      "pp.alignCounts",
      "pp.removeOutliers"
    ) %in%
      getNamespaceExports("playbase")
  ))

  set.seed(91)
  counts <- matrix(
    stats::rpois(30 * 8, 30),
    nrow = 30,
    dimnames = list(paste0("g", seq_len(30)), paste0("s", seq_len(8)))
  )
  counts[1, ] <- NA_real_
  counts[2, 3] <- NA_real_
  counts[3, 8] <- 1e7
  counts[4:30, 8] <- counts[4:30, 8] * seq(20, 200, length.out = 27)
  rownames(counts)[5:6] <- "dup"
  annot <- data.frame(
    feature = rownames(counts),
    row.names = seq_len(nrow(counts))
  )
  target <- stats::setNames(rep(c("A", "B"), each = 4), colnames(counts))
  batch <- data.frame(
    batch = rep(c("x", "y"), 4),
    row.names = colnames(counts)
  )
  samples <- data.frame(
    group = target,
    batch = batch$batch,
    row.names = colnames(counts)
  )
  contrasts <- matrix(
    target,
    ncol = 1,
    dimnames = list(colnames(counts), "A_vs_B")
  )
  options <- .opg_upload_preprocess_options(
    counts,
    datatype = "RNA-seq",
    filter_missing = TRUE,
    filter_threshold = 0.5,
    impute = TRUE,
    normalize = FALSE,
    remove_outliers = TRUE,
    outlier_threshold = 3,
    batch_correct = TRUE,
    batch_method = "limma",
    batch = batch,
    target = target,
    dedup = "average",
    max_features = 12L
  )

  first <- .opg_run_preprocess_preview(
    counts,
    samples,
    contrasts,
    annot,
    options,
    through = "final"
  )
  second <- .opg_run_preprocess_preview(
    first$counts,
    samples,
    contrasts,
    annot,
    first$options,
    through = "final"
  )

  testthat::expect_named(
    first,
    c(
      "counts",
      "X",
      "annot",
      "prior",
      "space",
      "alignment",
      "options"
    )
  )
  testthat::expect_identical(first$counts, counts)
  testthat::expect_identical(second$counts, counts)
  testthat::expect_identical(second$X, first$X)
  testthat::expect_identical(second$alignment, first$alignment)
  testthat::expect_lt(nrow(first$X), nrow(counts))
  testthat::expect_lt(ncol(first$X), ncol(counts))
  testthat::expect_identical(
    .opg_preview_counts(first),
    playbase::pp.alignCounts(first$counts, first$alignment, first$X)
  )
})

testthat::test_that("only current preprocessing metadata permits reanalysis", {
  current <- list(
    counts = matrix(1, 1, 1),
    settings = list(
      preprocess = list(
        options = list(normalize = TRUE),
        alignment = list(rows = list(1L), cols = 1L),
        space = "log2"
      )
    )
  )
  legacy <- list(
    counts = matrix(1, 1, 1),
    settings = list(norm_method = "CPM")
  )

  testthat::expect_true(.opg_has_preprocess_contract(current))
  testthat::expect_false(.opg_has_preprocess_contract(legacy))
})

testthat::test_that("reanalysis restores the canonical plural contrasts field", {
  contrasts <- matrix(
    c("A", "B"),
    ncol = 1,
    dimnames = list(c("s1", "s2"), "A_vs_B")
  )
  pgx <- list(
    organism = "human",
    samples = data.frame(group = c("A", "B"), row.names = c("s1", "s2")),
    counts = matrix(
      1:4,
      nrow = 2,
      dimnames = list(c("g1", "g2"), c("s1", "s2"))
    ),
    contrasts = contrasts,
    contrast = matrix("wrong", nrow = 2, ncol = 1),
    name = "rehydrated",
    description = "canonical upload payload",
    settings = list(
      preprocess = list(
        options = list(normalize = TRUE),
        alignment = list(rows = list(1L, 2L), cols = c(1L, 2L)),
        space = "log2"
      )
    )
  )

  testthat::expect_true(.opg_has_preprocess_contract(pgx))
  payload <- .opg_reanalysis_upload_payload(pgx)

  testthat::expect_identical(payload$contrasts.csv, contrasts)
  testthat::expect_identical(payload$samples.csv, pgx$samples)
  testthat::expect_identical(payload$counts.csv, pgx$counts)
  testthat::expect_identical(payload$organism, pgx$organism)
  testthat::expect_identical(payload$name, pgx$name)
  testthat::expect_identical(payload$description, pgx$description)
})
