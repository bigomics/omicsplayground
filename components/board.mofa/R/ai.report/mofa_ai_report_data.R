##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

MOFA_PROMPTS_DIR <- file.path(OPG, "components/board.mofa/prompts/mofa")

mofa_ai_compact_feature_table <- function(dt) {
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return("No top features available.")

  fmt <- data.frame(
    Feature = dt$feature,
    Symbol = dt$symbol,
    Weight = dt$weight,
    Centrality = dt$centrality,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  omicsai::omicsai_format_mdtable(
    fmt,
    formatters = list(
      Weight = function(x) omicsai::omicsai_format_num(x, 3),
      Centrality = function(x) omicsai::omicsai_format_num(x, 3)
    )
  )
}

mofa_ai_compact_pathway_table <- function(dt) {
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return("No significant pathway enrichment available.")

  cols <- intersect(c("pathway", "NES", "padj"), colnames(dt))
  dt <- dt[, cols, drop = FALSE]
  dt$pathway <- sub(".*:", "", dt$pathway)
  if ("padj" %in% colnames(dt)) colnames(dt)[colnames(dt) == "padj"] <- "q-value"

  omicsai::omicsai_format_mdtable(
    dt,
    formatters = list(
      NES = function(x) omicsai::omicsai_format_num(x, 3),
      `q-value` = omicsai::omicsai_format_pvalue
    )
  )
}

mofa_ai_build_summary_params <- function(mofa, pgx, factor_name, ntop = 12L) {
  ctx <- mofa_ai_extract_factor_data(mofa, pgx, factor_name, ntop = ntop)
  if (is.null(ctx)) return(NULL)

  dtype_summary <- "Not available"
  if (length(ctx$datatype_counts) > 0) {
    dtype_summary <- paste(sprintf("%s=%s", names(ctx$datatype_counts), as.integer(ctx$datatype_counts)), collapse = ", ")
  }

  list(
    factor = ctx$factor,
    traits = ctx$traits,
    experiment = ctx$experiment,
    headline_metrics = paste0(
      "Features: ", ctx$metrics$n_features,
      " | Significant pathways: ", ctx$metrics$n_sig_pathways,
      " | Max |NES|: ", omicsai::omicsai_format_num(ctx$metrics$max_abs_nes, 2),
      " | Max |weight|: ", omicsai::omicsai_format_num(ctx$metrics$max_abs_weight, 3)
    ),
    genesets = mofa_ai_compact_pathway_table(ctx$top_pathways),
    top_genes = mofa_ai_compact_feature_table(ctx$top_features),
    datatype_mix = dtype_summary
  )
}

mofa_ai_rank_factors <- function(factor_contexts) {
  if (length(factor_contexts) == 0) return(data.frame())

  rows <- lapply(factor_contexts, function(ctx) {
    score <-
      0.45 * min(ctx$metrics$n_sig_pathways / 10, 1) +
      0.30 * min(ctx$metrics$max_abs_nes / 2.5, 1) +
      0.25 * min(ctx$metrics$max_abs_weight / 1.5, 1)

    tier <- if (score >= 0.70) {
      "strong"
    } else if (score >= 0.45) {
      "moderate"
    } else {
      "weak"
    }

    data.frame(
      factor = ctx$factor,
      score = score,
      tier = tier,
      n_sig_pathways = ctx$metrics$n_sig_pathways,
      max_abs_nes = ctx$metrics$max_abs_nes,
      max_abs_weight = ctx$metrics$max_abs_weight,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out[order(-out$score), , drop = FALSE]
}

mofa_ai_build_report_tables <- function(mofa, pgx, max_contexts = 8L, ntop = 10L) {
  factors <- mofa_ai_factor_choices(mofa)
  if (length(factors) == 0) {
    return(list(text = "No MOFA factors available.", data = list(), ranking = data.frame()))
  }

  contexts <- list()
  for (k in factors) {
    ctx <- mofa_ai_extract_factor_data(mofa, pgx, k, ntop = ntop)
    if (!is.null(ctx)) contexts[[k]] <- ctx
  }

  if (length(contexts) == 0) {
    return(list(text = "No reportable MOFA factor contexts.", data = list(), ranking = data.frame()))
  }

  ranking <- mofa_ai_rank_factors(contexts)
  keep <- head(ranking$factor, as.integer(max_contexts))

  rank_table <- omicsai::omicsai_format_mdtable(
    head(ranking, as.integer(max_contexts)),
    formatters = list(
      score = function(x) omicsai::omicsai_format_num(x, 2),
      max_abs_nes = function(x) omicsai::omicsai_format_num(x, 2),
      max_abs_weight = function(x) omicsai::omicsai_format_num(x, 3)
    )
  )

  lines <- c(
    paste0("EXPERIMENT: ", mofa$experiment %||% pgx$name %||% pgx$description %||% "omics experiment"),
    "BOARD: MOFA",
    "",
    "## Factor Ranking",
    rank_table,
    "",
    "## Factor Detail"
  )

  for (k in keep) {
    ctx <- contexts[[k]]
    lines <- c(
      lines,
      paste0("### ", ctx$factor),
      paste0("- Traits: ", ctx$traits),
      paste0("- Significant pathways: ", ctx$metrics$n_sig_pathways,
        " | max |NES|: ", omicsai::omicsai_format_num(ctx$metrics$max_abs_nes, 2),
        " | max |weight|: ", omicsai::omicsai_format_num(ctx$metrics$max_abs_weight, 3)
      ),
      "",
      "Top pathways:",
      mofa_ai_compact_pathway_table(ctx$top_pathways),
      "",
      "Top weighted features:",
      mofa_ai_compact_feature_table(ctx$top_features),
      ""
    )
  }

  list(
    text = paste(lines, collapse = "\n"),
    data = contexts,
    ranking = ranking
  )
}

mofa_ai_build_methods <- function(mofa, pgx) {
  template <- omicsai::omicsai_load_template(
    file.path(MOFA_PROMPTS_DIR, "mofa_report_methods.md")
  )

  report_date <- format(Sys.Date(), "%Y-%m-%d")
  params <- list(
    experiment = mofa$experiment %||% pgx$name %||% pgx$description %||% "omics experiment",
    n_factors = length(mofa_ai_factor_choices(mofa)),
    n_samples = if (!is.null(mofa$F)) nrow(mofa$F) else NA_integer_,
    n_features = if (!is.null(mofa$W)) nrow(mofa$W) else NA_integer_,
    date = report_date
  )

  omicsai::collapse_lines(
    omicsai::omicsai_substitute_template(template, params),
    sprintf("_This report was generated with OmicsPlayground (BigOmics, %s)._", report_date),
    "_Note: AI-generated interpretation may contain inaccuracies and must be independently verified._",
    sep = "\n\n"
  )
}
