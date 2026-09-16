# Canonical preprocessing helpers for the bulk upload workflow.
#
# This file owns UI-to-engine option translation and staged previews.
# Numerical preprocessing always runs through the public leaf boundary.

# Extracts canonical row-layer labels from feature names.
# Matrices and data frames contribute their row names.
# Unprefixed features receive an empty label.
.opg_feature_layers <- function(x) {
  if (!is.null(dim(x))) {
    x <- rownames(x)
  }
  if (is.null(x)) {
    return(character())
  }
  ifelse(grepl(":", x), sub(":.*", "", x), "")
}

# Returns explicit layers only for a fully prefixed matrix.
# Prefixes must use the canonical alphanumeric `layer:feature` grammar.
# Single-omics matrices return NULL for whole-matrix dispatch.
.opg_preprocess_layers <- function(X) {
  feature_names <- rownames(X)
  if (
    is.null(feature_names) ||
      !length(feature_names) ||
      !all(grepl("^[A-Za-z0-9]+:", feature_names))
  ) {
    return(NULL)
  }
  sub(":.*", "", feature_names)
}

# Imputes a display or analysis matrix through the canonical family API.
# Fully prefixed matrices are processed independently by layer.
# The function performs no pgx-object or option translation.
.opg_impute <- function(X, method = "SVD2") {
  playbase.preprocess::pp.impute(
    as.matrix(X),
    layers = .opg_preprocess_layers(X),
    method = method
  )
}

# Converts a count matrix to canonical log-CPM values.
# Conversion and normalization both use the declared additive prior.
# Sparse inputs are materialized because the leaf API is matrix-only.
.opg_log_cpm <- function(counts, total = 1e6, prior = 1) {
  counts <- as.matrix(counts)
  X <- playbase.preprocess::pp.convertSpace(
    counts,
    from = "counts",
    to = "log2",
    prior = prior
  )
  X <- playbase.preprocess::pp.normalize(
    X,
    method = "CPM",
    space = "log2",
    prior = prior
  )
  if (!identical(total, 1e6)) {
    scaled <- playbase.preprocess::pp.convertSpace(
      X,
      from = "log2",
      to = "counts",
      prior = prior
    ) *
      (total / 1e6)
    X <- playbase.preprocess::pp.convertSpace(
      scaled,
      from = "counts",
      to = "log2",
      prior = prior
    )
  }
  X
}

# Normalizes linear counts for the UI comparison panel.
# Only canonical matrix-family methods are offered by the panel.
# The output is returned to linear count scale for its existing plots.
.opg_count_normalize_preview <- function(counts, method = "CPM") {
  counts <- as.matrix(counts)
  if (identical(method, "none")) {
    return(counts)
  }
  X <- playbase.preprocess::pp.convertSpace(
    counts,
    from = "counts",
    to = "log2",
    prior = 1
  )
  normalize_args <- list(
    X = X,
    method = method,
    space = "log2"
  )
  if (method %in% c("CPM", "CPM+quantile", "TMM")) {
    normalize_args$prior <- 1
  }
  X <- do.call(playbase.preprocess::pp.normalize, normalize_args)
  playbase.preprocess::pp.convertSpace(
    X,
    from = "log2",
    to = "counts",
    prior = 1
  )
}

# Builds the canonical options for one bulk upload.
# Datatype policy is translated once before calling the matrix engine.
# Batch inputs remain plain sample-aligned data in the options list.
.opg_upload_preprocess_options <- function(
  counts,
  datatype,
  is_npx = FALSE,
  zero_as_na = FALSE,
  filter_missing = FALSE,
  filter_threshold = 3,
  impute = FALSE,
  impute_method = "SVD2",
  normalize = TRUE,
  norm_method = "CPM",
  ref_gene = NULL,
  remove_outliers = FALSE,
  outlier_threshold = 3,
  batch_correct = FALSE,
  batch_method = "limma",
  batch = NULL,
  target = NULL,
  meth_type = NULL,
  dedup = NULL,
  max_features = NULL
) {
  layers <- .opg_preprocess_layers(counts)
  input_space <- "counts"
  output_space <- "log2"
  normalize_args <- list()

  if (identical(datatype, "methylomics")) {
    input_space <- "beta"
    output_space <- "beta"
    impute <- FALSE
    if (identical(norm_method, "BMIQ")) {
      normalize_args <- list(meth_type = meth_type)
    }
  } else if (!is.null(layers)) {
    present <- unique(layers)
    norm_method <- stats::setNames(
      ifelse(present == "gx", "CPM", "maxMedian"),
      present
    )
  } else if (identical(norm_method, "reference")) {
    normalize_args <- list(ref = ref_gene)
  }

  options <- list(
    input_space = input_space,
    output_space = output_space,
    layers = layers,
    is_npx = isTRUE(is_npx),
    zero_as_na = isTRUE(zero_as_na),
    filter_missing = isTRUE(filter_missing),
    filter_threshold = as.numeric(filter_threshold),
    impute = isTRUE(impute),
    impute_method = impute_method,
    impute_args = list(),
    normalize = isTRUE(normalize),
    norm_method = norm_method,
    normalize_args = normalize_args,
    remove_outliers = isTRUE(remove_outliers),
    outlier_threshold = as.numeric(outlier_threshold),
    outlier_methods = c("z.correlation", "z.distance", "z.features"),
    batch_correct = isTRUE(batch_correct),
    batch_method = batch_method,
    batch = batch,
    target = target,
    batch_args = list(),
    max_features = max_features
  )
  if (!is.null(dedup)) {
    options$dedup <- dedup
  }
  options
}

# Derives sample-aligned target and selected batch columns for the upload UI.
# Autodetection remains playbase policy and explicit selections retain order.
# Disabled correction returns NULL inputs without running batch diagnostics.
.opg_upload_batch_inputs <- function(
  counts,
  samples,
  contrasts,
  selection,
  enabled
) {
  if (!isTRUE(enabled)) {
    return(list(target = NULL, batch = NULL))
  }
  sample_names <- colnames(counts)
  if (!is.null(sample_names) && !is.null(rownames(samples))) {
    if (!all(sample_names %in% rownames(samples))) {
      stop("Uploaded samples do not cover every count-matrix column")
    }
    samples <- samples[sample_names, , drop = FALSE]
  }
  if (!is.null(sample_names) && !is.null(rownames(contrasts))) {
    if (!all(sample_names %in% rownames(contrasts))) {
      stop("Uploaded contrasts do not cover every count-matrix column")
    }
    contrasts <- contrasts[sample_names, , drop = FALSE]
  }

  target <- playbase::contrasts2pheno(contrasts, samples)
  selected <- selection
  if (is.null(selected) || any(selected %in% "<none>")) {
    selected <- character()
  } else if (any(selected %in% "<autodetect>")) {
    selected <- playbase::get_model_parameters(
      as.matrix(counts),
      samples,
      pheno = NULL,
      contrasts = contrasts
    )$batch.pars
  }
  selected <- intersect(selected, colnames(samples))
  batch <- if (length(selected)) {
    samples[, selected, drop = FALSE]
  } else {
    NULL
  }
  list(target = target, batch = batch)
}

# Disables every pipeline family after a requested preview stage.
# The retained prefix always follows the production hardcoded order.
# Final previews keep the caller's complete canonical options unchanged.
.opg_preview_options <- function(options, through = "final") {
  stages <- c("impute", "normalize", "outliers", "batch", "final")
  through <- match.arg(through, stages)
  if (match(through, stages) < match("normalize", stages)) {
    options$normalize <- FALSE
  }
  if (match(through, stages) < match("outliers", stages)) {
    options$remove_outliers <- FALSE
  }
  if (match(through, stages) < match("batch", stages)) {
    options$batch_correct <- FALSE
  }
  if (match(through, stages) < match("final", stages)) {
    options$dedup <- "off"
    options$max_features <- NULL
  }
  options
}

# Runs one staged upload preview through the real preprocessing boundary.
# Samples and contrasts enter the leaf as plain aligned metadata.
# The returned value is the canonical seven-field preprocessing result.
.opg_run_preprocess_preview <- function(
  counts,
  samples,
  contrasts,
  annot,
  options,
  through = "final"
) {
  playbase.preprocess::pgx.preprocess(
    counts = as.matrix(counts),
    samples = samples,
    contrasts = contrasts,
    annot = annot,
    options = .opg_preview_options(options, through = through)
  )
}

# Aligns pristine source counts to one canonical preview result.
# Positional alignment handles filtering, deduplication, and outlier columns.
# The returned matrix always has the preview X dimensions and names.
.opg_preview_counts <- function(result) {
  playbase.preprocess::pp.alignCounts(
    result$counts,
    result$alignment,
    X = result$X
  )
}

# Aligns an earlier processed preview to a later preview's final positions.
# Row groups and source-column indices are matched without feature-name joins.
# Every requested target position must exist in the earlier preview.
.opg_align_preview_X <- function(result, target) {
  row_key <- vapply(result$alignment$rows, paste, collapse = ",", character(1))
  target_row_key <- vapply(
    target$alignment$rows,
    paste,
    collapse = ",",
    character(1)
  )
  rows <- match(target_row_key, row_key)
  cols <- match(target$alignment$cols, result$alignment$cols)
  if (anyNA(rows) || anyNA(cols)) {
    stop("Preview stages do not share a valid positional alignment")
  }
  out <- result$X[rows, cols, drop = FALSE]
  dimnames(out) <- dimnames(target$X)
  out
}

# Checks whether a pgx object can safely seed canonical reanalysis.
# Current objects carry resolved options and direct alignment metadata.
# Legacy normalization records are deliberately not interpreted.
.opg_has_preprocess_contract <- function(pgx) {
  state <- pgx$settings$preprocess
  is.list(state) &&
    is.list(state$options) &&
    is.list(state$alignment) &&
    !is.null(state$space) &&
    !is.null(pgx$counts)
}

# Returns the canonical upload values restored during pgx reanalysis.
# Contrast data always comes from the plural pgx `contrasts` field.
# Callers must verify the preprocessing contract before using the payload.
.opg_reanalysis_upload_payload <- function(pgx) {
  list(
    organism = pgx$organism,
    samples.csv = pgx$samples,
    contrasts.csv = pgx$contrasts,
    counts.csv = pgx$counts,
    name = pgx$name,
    description = pgx$description
  )
}

# Aligns a pgx object's pristine source matrix to its analysis matrix.
# Current preprocessing metadata is mandatory; names are never a fallback.
# The independent single-cell path retains its explicit name-aligned contract.
.opg_pgx_analysis_counts <- function(pgx) {
  if (!is.null(pgx$datatype) && pgx$datatype %in% c("scRNA-seq", "scRNAseq")) {
    return(pgx$counts[rownames(pgx$X), colnames(pgx$X), drop = FALSE])
  }
  if (!.opg_has_preprocess_contract(pgx)) {
    stop("Dataset has no canonical preprocessing alignment")
  }
  playbase.preprocess::pp.alignCounts(
    pgx$counts,
    pgx$settings$preprocess$alignment,
    X = pgx$X
  )
}

# Applies the upload-only nonzero-median preview normalization.
# This method is not part of the production preprocessing vocabulary.
# It remains local because it only conditions a sample-preview graphic.
.opg_median_center_nonzero <- function(X) {
  X <- as.matrix(X)
  centers <- apply(X, 2L, function(x) stats::median(x[x > 0], na.rm = TRUE))
  t(t(X) / (1e-8 + centers)) * stats::median(centers, na.rm = TRUE)
}
