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

# Imputes a display or analysis matrix through the canonical family API.
# Fully prefixed matrices are processed independently by layer.
# The function performs no pgx-object or option translation.
.opg_impute <- function(X, method = "SVD2") {
  playbase.preprocess::pp.impute(
    as.matrix(X),
    layers = playbase.preprocess::pp.inferLayers(X),
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
  batch.correct.method = "no_batch_correct",
  batch = NULL,
  target = NULL,
  meth_type = NULL,
  dedup = NULL,
  max_features = NULL
) {
  layers <- playbase.preprocess::pp.inferLayers(counts)
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

  ## ComBat and limma need an explicit batch variable. Autodetection can
  ## resolve nothing, so the stage stays disabled instead of erroring.
  if (
    any(batch.correct.method %in% c("ComBat", "limma")) &&
      (is.null(batch) || !NCOL(batch))
  ) {
    batch.correct.method <- "no_batch_correct"
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
    batch.correct.method = batch.correct.method,
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
    options$batch.correct.method <- "no_batch_correct"
  }
  if (match(through, stages) < match("final", stages)) {
    options$dedup <- "off"
    options$max_features <- NULL
  }
  options
}

# Splits canonical upload options into the two pgx.createPGX batch inputs.
# Playbase owns the batch selector as a top-level argument, never a nested key.
# Resolved batch columns are reported so Playbase never re-autodetects them.
.opg_createpgx_preprocess <- function(options) {
  if (is.null(options)) {
    return(list(
      preprocess = NULL,
      batch.correct.method = "no_batch_correct",
      batch.pars = "<none>"
    ))
  }
  method <- options$batch.correct.method
  if (is.null(method)) method <- "no_batch_correct"
  options$batch.correct.method <- NULL
  batch.pars <- colnames(options$batch)
  if (!length(batch.pars)) batch.pars <- "<none>"
  list(
    preprocess = options,
    batch.correct.method = unname(method[[1L]]),
    batch.pars = batch.pars
  )
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

# Runs one staged preview and reports engine failures to the open panel.
# Method preconditions such as confounded batches are user-facing data states.
# A successful run returns the canonical preprocessing result unchanged.
.opg_validated_preview <- function(
  counts,
  samples,
  contrasts,
  annot,
  options,
  through = "final",
  label = "Preprocessing"
) {
  failure <- NULL
  result <- tryCatch(
    .opg_run_preprocess_preview(
      counts = counts,
      samples = samples,
      contrasts = contrasts,
      annot = annot,
      options = options,
      through = through
    ),
    error = function(e) {
      failure <<- conditionMessage(e)
      NULL
    }
  )
  if (!is.null(failure)) {
    shiny::validate(shiny::need(FALSE, paste0(label, " failed: ", failure)))
  }
  result
}

# Averages duplicate rows in a canonical preview result.
# Declared result spaces and inferred layers select each averaging rule.
# The returned matrix preserves every non-duplicate preview row.
.opg_deduplicate_preview <- function(result) {
  X <- result$X
  if (!anyDuplicated(rownames(X))) {
    return(X)
  }
  playbase.preprocess::pp.deduplicate(
    X,
    method = "average",
    space = result$space,
    layers = playbase.preprocess::pp.inferLayers(X)
  )$X
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
# Current objects use positional metadata, while single-cell objects use exact
# names. Historical bulk objects without alignment are rejected explicitly.
.opg_pgx_analysis_counts <- function(pgx) {
  state <- pgx$settings$preprocess
  if (is.list(state) && is.list(state$alignment)) {
    return(playbase.preprocess::pp.alignCounts(
      pgx$counts,
      state$alignment,
      X = pgx$X
    ))
  }
  single_cell <- !is.null(pgx$datatype) &&
    pgx$datatype %in% c("scRNA-seq", "scRNAseq")

  # Datasets built before the canonical alignment, and single-cell datasets that
  # never carry one, are positioned by matching axis names instead.
  aligned <- .opg_name_aligned_counts(
    pgx$counts,
    pgx$X,
    strict = single_cell
  )
  if (!is.null(aligned)) {
    return(aligned)
  }
  # A legacy dataset whose axes cannot be matched keeps its source counts, which
  # is what every board consumed before the alignment existed.
  pgx$counts
}

# Positions source counts at the analysis axes using row and column names.
# Strict callers require a complete, unique, and exact name correspondence.
# Lenient callers receive NULL when the axes cannot be matched by name.
.opg_name_aligned_counts <- function(counts, X, strict) {
  refuse <- function(message) {
    if (strict) {
      stop(message, call. = FALSE)
    }
    NULL
  }
  if (is.null(dim(counts)) || is.null(dim(X))) {
    return(refuse("Single-cell dataset counts and X must both be matrices"))
  }
  source_rows <- rownames(counts)
  source_cols <- colnames(counts)
  analysis_rows <- rownames(X)
  analysis_cols <- colnames(X)
  axes <- list(source_rows, source_cols, analysis_rows, analysis_cols)
  valid_axes <- vapply(
    axes,
    function(x) {
      !is.null(x) && !anyNA(x) && all(nzchar(x))
    },
    logical(1)
  )
  if (!all(valid_axes)) {
    return(refuse(
      "Single-cell dataset counts and X must have complete row and column names"
    ))
  }
  if (
    anyDuplicated(source_rows) ||
      anyDuplicated(source_cols) ||
      anyDuplicated(analysis_rows) ||
      anyDuplicated(analysis_cols)
  ) {
    return(refuse(
      "Single-cell dataset counts and X must have unique axis names"
    ))
  }
  rows <- match(analysis_rows, source_rows)
  cols <- match(analysis_cols, source_cols)
  if (anyNA(rows) || anyNA(cols)) {
    return(refuse(
      "Single-cell dataset X row and column names must exist exactly in counts"
    ))
  }
  counts[rows, cols, drop = FALSE]
}

# Applies the upload-only nonzero-median preview normalization.
# This method is not part of the production preprocessing vocabulary.
# It remains local because it only conditions a sample-preview graphic.
.opg_median_center_nonzero <- function(X) {
  X <- as.matrix(X)
  centers <- apply(X, 2L, function(x) stats::median(x[x > 0], na.rm = TRUE))
  t(t(X) / (1e-8 + centers)) * stats::median(centers, na.rm = TRUE)
}
