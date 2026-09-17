# Canonical bulk-upload preprocessing contract tests.
#
# These tests cover staged previews, pristine source reuse, and alignment.
# They exercise the installed preprocessing leaf used by the application.

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
      "batch.correct.method",
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
  testthat::expect_identical(
    .opg_preview_options(options, "outliers")$batch.correct.method,
    "no_batch_correct"
  )
})

testthat::test_that("upload layer inference preserves sparse inputs", {
  counts <- Matrix::sparseMatrix(
    i = c(1L, 2L, 3L, 4L),
    j = c(1L, 2L, 1L, 2L),
    x = c(1, 2, 3, 4),
    dims = c(4L, 2L),
    dimnames = list(
      c("gx:a", "gx:b", "px:a", "px:b"),
      c("s1", "s2")
    )
  )
  observed_class <- NULL
  infer_layers <- playbase.preprocess::pp.inferLayers
  testthat::local_mocked_bindings(
    pp.inferLayers = function(X) {
      observed_class <<- class(X)
      infer_layers(X)
    },
    .package = "playbase.preprocess"
  )

  options <- .opg_upload_preprocess_options(
    counts,
    datatype = "multi-omics",
    normalize = FALSE
  )

  testthat::expect_s4_class(counts, "dgCMatrix")
  testthat::expect_true("dgCMatrix" %in% observed_class)
  testthat::expect_identical(options$layers, c("gx", "gx", "px", "px"))
  testthat::expect_identical(
    options$norm_method,
    c(gx = "CPM", px = "maxMedian")
  )
})

testthat::test_that("duplicate previews average in declared log2 space", {
  X <- matrix(
    c(0, 2, 2, 4, 1, 3),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("dup", "dup", "solo"), c("s1", "s2"))
  )

  result <- .opg_deduplicate_preview(list(X = X, space = "log2"))

  testthat::expect_identical(rownames(result), c("dup", "solo"))
  testthat::expect_equal(
    unname(result["dup", ]),
    unname(log2(colMeans(2^X[1:2, , drop = FALSE])))
  )
  testthat::expect_identical(result["solo", ], X["solo", ])
})

testthat::test_that("duplicate previews honor beta and layer spaces", {
  X <- matrix(
    c(0.2, 0.8, 0.8, 0.2, 0, 2, 2, 4),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      c("mx:dup", "mx:dup", "gx:dup", "gx:dup"),
      c("s1", "s2")
    )
  )

  result <- .opg_deduplicate_preview(list(
    X = X,
    space = c(mx = "beta", gx = "log2")
  ))

  testthat::expect_identical(rownames(result), c("mx:dup", "gx:dup"))
  testthat::expect_equal(unname(result["mx:dup", ]), c(0.5, 0.5))
  testthat::expect_equal(
    unname(result["gx:dup", ]),
    unname(log2(colMeans(2^X[3:4, , drop = FALSE])))
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
      getNamespaceExports("playbase.preprocess")
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
    batch.correct.method = "limma",
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
  direct <- playbase.preprocess::pgx.preprocess(
    counts = counts,
    samples = samples,
    contrasts = contrasts,
    annot = annot,
    options = .opg_preview_options(options, through = "final")
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
  testthat::expect_identical(first, direct)
  testthat::expect_identical(second$counts, counts)
  testthat::expect_identical(second$X, first$X)
  testthat::expect_identical(second$alignment, first$alignment)
  testthat::expect_lt(nrow(first$X), nrow(counts))
  testthat::expect_lt(ncol(first$X), ncol(counts))
  testthat::expect_identical(
    .opg_preview_counts(first),
    playbase.preprocess::pp.alignCounts(
      first$counts,
      first$alignment,
      first$X
    )
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

testthat::test_that("analysis counts use canonical positional metadata", {
  counts <- matrix(
    seq_len(9),
    nrow = 3,
    dimnames = list(c("g1", "g2", "g3"), c("s1", "s2", "s3"))
  )
  X <- matrix(
    0,
    nrow = 2,
    ncol = 2,
    dimnames = list(
      c("processed-g2", "processed-g1"),
      c("sample-3", "sample-1")
    )
  )
  pgx <- list(
    counts = counts,
    X = X,
    settings = list(
      preprocess = list(
        options = list(normalize = TRUE),
        alignment = list(rows = list(2L, 1L), cols = c(3L, 1L)),
        space = "log2"
      )
    )
  )
  expected <- counts[c(2L, 1L), c(3L, 1L), drop = FALSE]
  dimnames(expected) <- dimnames(X)

  testthat::expect_identical(.opg_pgx_analysis_counts(pgx), expected)

  pgx$settings$preprocess$options <- NULL
  testthat::expect_false(.opg_has_preprocess_contract(pgx))
  testthat::expect_identical(.opg_pgx_analysis_counts(pgx), expected)
})

testthat::test_that("analysis counts align legacy and single-cell names", {
  counts <- matrix(
    seq_len(12),
    nrow = 3,
    dimnames = list(c("g1", "g2", "g3"), c("s1", "s2", "s3", "s4"))
  )
  X <- matrix(
    0,
    nrow = 2,
    ncol = 2,
    dimnames = list(c("g3", "g1"), c("s4", "s2"))
  )
  legacy <- list(
    counts = counts,
    X = X,
    settings = list(norm_method = "CPM")
  )
  expected <- counts[c("g3", "g1"), c("s4", "s2"), drop = FALSE]

  ## Datasets built before the alignment existed are positioned by name.
  testthat::expect_identical(.opg_pgx_analysis_counts(legacy), expected)

  for (datatype in c("scRNA-seq", "scRNAseq")) {
    legacy$datatype <- datatype
    testthat::expect_identical(.opg_pgx_analysis_counts(legacy), expected)
  }

  ## Single-cell datasets still require an exact name correspondence.
  unmatched <- legacy
  rownames(unmatched$X)[1L] <- "missing"
  testthat::expect_error(
    .opg_pgx_analysis_counts(unmatched),
    "must exist exactly in counts"
  )

  ## A legacy bulk dataset that cannot be matched keeps its source counts.
  unmatched$datatype <- "RNA-seq"
  testthat::expect_identical(.opg_pgx_analysis_counts(unmatched), counts)

  unnamed <- list(counts = unname(counts), X = unname(X), datatype = "RNA-seq")
  testthat::expect_identical(
    .opg_pgx_analysis_counts(unnamed),
    unname(counts)
  )
})

testthat::test_that("background creation uses the canonical preprocess handoff", {
  process_file <- file.path(.opg_test_root, "bin", "pgxcreate_op.R")
  process_source <- paste(readLines(process_file), collapse = "\n")

  testthat::expect_true(grepl(
    "preprocess = params$preprocess",
    process_source,
    fixed = TRUE
  ))
  ## Batch selection is a top-level Playbase argument, never a nested key.
  testthat::expect_true(grepl(
    "batch.correct.method = batch.correct.method",
    process_source,
    fixed = TRUE
  ))
  testthat::expect_true(grepl(
    "batch.pars = batch.pars",
    process_source,
    fixed = TRUE
  ))
})

testthat::test_that("supervised correction without batch stays disabled", {
  counts <- matrix(
    seq_len(8),
    nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("s", seq_len(4)))
  )
  batch <- data.frame(
    plate = rep(c("x", "y"), 2),
    row.names = colnames(counts)
  )

  for (method in c("ComBat", "limma")) {
    options <- .opg_upload_preprocess_options(
      counts,
      datatype = "RNA-seq",
      batch.correct.method = method
    )
    testthat::expect_identical(
      options$batch.correct.method,
      "no_batch_correct"
    )
    kept <- .opg_upload_preprocess_options(
      counts,
      datatype = "RNA-seq",
      batch.correct.method = method,
      batch = batch
    )
    testthat::expect_identical(kept$batch.correct.method, method)
  }

  ## Latent-factor methods estimate their own structure and stay selected.
  for (method in c("RUV", "SVA", "NPM")) {
    options <- .opg_upload_preprocess_options(
      counts,
      datatype = "RNA-seq",
      batch.correct.method = method
    )
    testthat::expect_identical(options$batch.correct.method, method)
  }
})

testthat::test_that("preview failures surface as panel validation", {
  counts <- matrix(
    rep(seq_len(6), 4) + rep(c(0, 1, 5, 6), each = 6),
    nrow = 6,
    dimnames = list(paste0("g", seq_len(6)), paste0("s", seq_len(4)))
  )
  target <- stats::setNames(rep(c("A", "B"), each = 2), colnames(counts))
  ## Batch tracks the target exactly, which ComBat cannot adjust.
  batch <- data.frame(
    site = rep(c("x", "y"), each = 2),
    row.names = colnames(counts)
  )
  samples <- data.frame(
    group = target,
    site = batch$site,
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
    batch.correct.method = "ComBat",
    batch = batch,
    target = target
  )

  ## A batch confounded with the target falls back to the unprotected fit,
  ## which is the production behaviour the upload flow relies on.
  testthat::expect_message(
    corrected <- .opg_run_preprocess_preview(
      counts,
      samples,
      contrasts,
      annot = NULL,
      options = options,
      through = "batch"
    ),
    "without target protection"
  )
  testthat::expect_true(all(is.finite(corrected$X)))

  ## A method precondition that cannot be satisfied reaches the panel.
  broken <- .opg_upload_preprocess_options(
    counts,
    datatype = "RNA-seq",
    batch.correct.method = "SVA",
    target = stats::setNames(rep("A", ncol(counts)), colnames(counts))
  )
  testthat::expect_error(
    .opg_run_preprocess_preview(
      counts,
      samples,
      contrasts,
      annot = NULL,
      options = broken,
      through = "batch"
    ),
    "at least two groups"
  )
  testthat::expect_error(
    .opg_validated_preview(
      counts,
      samples,
      contrasts,
      annot = NULL,
      options = broken,
      through = "batch",
      label = "Batch-effect correction"
    ),
    class = "shiny.silent.error"
  )
})

testthat::test_that("createPGX handoff splits out the batch selector", {
  counts <- matrix(
    seq_len(8),
    nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("s", seq_len(4)))
  )
  batch <- data.frame(
    plate = rep(c("x", "y"), 2),
    row.names = colnames(counts)
  )
  options <- .opg_upload_preprocess_options(
    counts,
    datatype = "RNA-seq",
    batch.correct.method = "limma",
    batch = batch
  )

  handoff <- .opg_createpgx_preprocess(options)
  testthat::expect_identical(handoff$batch.correct.method, "limma")
  testthat::expect_identical(handoff$batch.pars, "plate")
  testthat::expect_false("batch.correct.method" %in% names(handoff$preprocess))
  testthat::expect_identical(handoff$preprocess$batch, batch)

  disabled <- .opg_createpgx_preprocess(
    .opg_upload_preprocess_options(counts, datatype = "RNA-seq")
  )
  testthat::expect_identical(
    disabled$batch.correct.method,
    "no_batch_correct"
  )
  testthat::expect_identical(disabled$batch.pars, "<none>")

  absent <- .opg_createpgx_preprocess(NULL)
  testthat::expect_null(absent$preprocess)
  testthat::expect_identical(absent$batch.correct.method, "no_batch_correct")
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
