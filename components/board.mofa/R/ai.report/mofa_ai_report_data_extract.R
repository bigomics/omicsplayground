##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

mofa_ai_factor_choices <- function(mofa) {
  if (is.null(mofa$W) || is.null(colnames(mofa$W))) return(character(0))
  colnames(mofa$W)
}

mofa_ai_trait_summary <- function(mofa, factor_name) {
  if (is.null(mofa$Y) || is.null(mofa$F) || !factor_name %in% colnames(mofa$F)) {
    return("None")
  }

  factor_scores <- mofa$F[, factor_name]
  trait_names <- colnames(mofa$Y)
  if (length(trait_names) == 0) return("None")

  cors <- sapply(trait_names, function(t) {
    y <- mofa$Y[, t]
    if (!is.numeric(y) || sum(!is.na(y)) < 4) return(NA_real_)
    tryCatch(cor(factor_scores, y, use = "pairwise.complete.obs"), error = function(e) NA_real_)
  })
  cors <- cors[!is.na(cors)]

  if (length(cors) == 0) return("None")

  sig <- cors[abs(cors) > 0.3]
  picked <- if (length(sig) > 0) {
    sig[order(-abs(sig))]
  } else {
    cors[order(-abs(cors))][seq_len(min(3L, length(cors)))]
  }

  paste(
    sprintf("%s (r=%s)", names(picked), omicsai::omicsai_format_num(picked, 2)),
    collapse = ", "
  )
}

mofa_ai_extract_factor_data <- function(mofa, pgx, factor_name, ntop = 12L) {
  if (is.null(mofa) || is.null(factor_name) || !nzchar(factor_name)) return(NULL)
  if (is.null(mofa$W) || !factor_name %in% colnames(mofa$W)) return(NULL)

  w <- mofa$W[, factor_name]
  w <- w[!is.na(w)]
  if (length(w) == 0) return(NULL)

  top_idx <- multiomics_ai_top_idx_by_abs(ntop, w)
  top_features <- names(w)[top_idx]
  top_weights <- w[top_features]

  symbols <- top_features
  if (!is.null(pgx$genes)) {
    symbols <- tryCatch(
      playbase::probe2symbol(top_features, pgx$genes, "symbol", fill_na = TRUE),
      error = function(e) top_features
    )
    symbols <- ifelse(is.na(symbols) | symbols == "", top_features, symbols)
  }

  centrality_vals <- rep(NA_real_, length(top_features))
  centrality_data <- tryCatch(
    playbase::mofa.plot_centrality(mofa, factor_name, justdata = TRUE),
    error = function(e) NULL
  )
  if (is.data.frame(centrality_data) &&
    "centrality" %in% colnames(centrality_data) &&
    "feature" %in% colnames(centrality_data)) {
    idx <- match(top_features, centrality_data$feature)
    centrality_vals <- centrality_data$centrality[idx]
  }

  prefix <- tryCatch(playbase::mofa.get_prefix(top_features), error = function(e) NULL)
  dtype_counts <- if (!is.null(prefix)) sort(table(prefix), decreasing = TRUE) else table(character(0))

  pathways <- data.frame()
  if (!is.null(mofa$gsea) && !is.null(mofa$gsea$table) && factor_name %in% names(mofa$gsea$table)) {
    g <- mofa$gsea$table[[factor_name]]
    if (is.data.frame(g) && nrow(g) > 0 && all(c("pathway", "NES", "padj") %in% colnames(g))) {
      g <- g[order(-abs(g$NES)), , drop = FALSE]
      pathways <- g[!is.na(g$padj) & g$padj < 0.05, , drop = FALSE]
      pathways <- head(pathways, 10L)
    }
  }

  max_abs_nes <- if (nrow(pathways) > 0) max(abs(pathways$NES), na.rm = TRUE) else 0
  max_abs_weight <- max(abs(top_weights), na.rm = TRUE)

  list(
    factor = factor_name,
    experiment = multiomics_ai_experiment_label(pgx, mofa$experiment),
    traits = mofa_ai_trait_summary(mofa, factor_name),
    top_features = data.frame(
      feature = top_features,
      symbol = symbols,
      weight = as.numeric(top_weights),
      centrality = as.numeric(centrality_vals),
      stringsAsFactors = FALSE
    ),
    top_pathways = pathways,
    datatype_counts = dtype_counts,
    metrics = list(
      n_features = length(w),
      n_sig_pathways = nrow(pathways),
      max_abs_nes = max_abs_nes,
      max_abs_weight = max_abs_weight
    )
  )
}
