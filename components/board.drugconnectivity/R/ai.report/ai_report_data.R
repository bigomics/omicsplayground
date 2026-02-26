##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

DRUGCONNECTIVITY_PROMPTS_DIR <- file.path(OPG, "components/board.drugconnectivity/prompts")

dc_percent <- function(x, digits = 1L) {
  sprintf(paste0("%.", digits, "f%%"), 100 * x)
}

dc_analysis_type_info <- function(analysis_type) {
  desc <- switch(
    analysis_type,
    "L1000_ACTIVITYS_N20D1011" = paste(
      "L1000 activity-signature connectivity analysis.",
      "Experimental signatures are matched against perturbational transcriptional profiles",
      "from LINCS L1000 to quantify mimicking vs opposing drug connectivity."
    ),
    "CTRPv2/Sensitivity" = paste(
      "CTRPv2 sensitivity-based drug connectivity analysis.",
      "Drug associations are interpreted in the context of pharmacologic sensitivity",
      "profiles from the Cancer Therapeutics Response Portal v2."
    ),
    "GDSC/Sensitivity" = paste(
      "GDSC sensitivity-based drug connectivity analysis.",
      "Drug associations are interpreted in the context of pharmacogenomic sensitivity",
      "profiles from the Genomics of Drug Sensitivity in Cancer resource."
    ),
    paste(
      "Drug connectivity analysis based on the selected pre-computed signature resource.",
      "Interpret results according to the selected analysis type metadata."
    )
  )

  list(analysis_type = analysis_type %||% "unknown", analysis_type_description = desc)
}

dc_compact_drug_table <- function(dt, n = 10L) {
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return("No data available.")
  dt <- head(dt, n)
  fmt <- data.frame(
    Drug = dt$drug,
    NES = dt$NES,
    `q-value` = dt$padj,
    MOA = ifelse(is.na(dt$moa) | dt$moa == "", "-", dt$moa),
    Target = ifelse(is.na(dt$target) | dt$target == "", "-", dt$target),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  omicsai::omicsai_format_mdtable(
    fmt,
    formatters = list(
      NES = function(x) omicsai::omicsai_format_num(x, 3),
      `q-value` = omicsai::omicsai_format_pvalue
    )
  )
}

dc_compact_moa_table <- function(dt, n = 10L, label = "Pathway") {
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return("No significant enrichment available.")
  dt <- head(dt, n)
  fmt <- data.frame(
    Name = dt$pathway,
    NES = dt$NES,
    `q-value` = dt$padj,
    Size = dt$size,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(fmt)[1] <- label
  omicsai::omicsai_format_mdtable(
    fmt,
    formatters = list(
      NES = function(x) omicsai::omicsai_format_num(x, 3),
      `q-value` = omicsai::omicsai_format_pvalue
    )
  )
}

drugconnectivity_build_summary_params <- function(pgx,
                                                  contrast,
                                                  method,
                                                  only_annotated = FALSE,
                                                  ntop = 12L) {
  ctx <- extract_drugconnectivity_context_data(
    pgx = pgx,
    contrast = contrast,
    method = method,
    only_annotated = only_annotated,
    n_top = ntop
  )
  if (is.null(ctx)) return(NULL)

  rel <- ctx$reliability
  at <- dc_analysis_type_info(method)
  headline_metrics <- paste0(
    "Drugs tested: ", ctx$n_drugs_tested,
    " | significant (q<0.05): ", ctx$n_drugs_sig,
    " | opposing: ", ctx$n_neg,
    " | mimicking: ", ctx$n_pos,
    " | max |NES|: ", omicsai::omicsai_format_num(ctx$max_abs_NES, 2)
  )

  reliability_notes <- paste(
    paste0("Annotated drug coverage: ", dc_percent(rel$frac_annotated)),
    paste0("Significant MOA classes: ", rel$n_sig_moa_classes),
    paste0("Significant targets: ", rel$n_sig_targets),
    paste(ctx$caveats, collapse = " "),
    sep = "  \\n"
  )

  list(
    contrast = ctx$contrast,
    method = ctx$method,
    analysis_type = at$analysis_type,
    analysis_type_description = at$analysis_type_description,
    experiment = ctx$experiment,
    headline_metrics = headline_metrics,
    top_opposing_table = dc_compact_drug_table(ctx$top_opposing, n = ntop),
    top_mimicking_table = dc_compact_drug_table(ctx$top_mimicking, n = ntop),
    moa_class_table = dc_compact_moa_table(ctx$moa_class, n = 10L, label = "MOA Class"),
    moa_target_table = dc_compact_moa_table(ctx$moa_target, n = 10L, label = "Target"),
    reliability_notes = reliability_notes
  )
}

drugconnectivity_rank_contexts <- function(context_data) {
  if (length(context_data) == 0) return(data.frame())

  rows <- lapply(context_data, function(ctx) {
    rel <- ctx$reliability
    score <-
      0.35 * min(ctx$n_drugs_sig / 25, 1) +
      0.25 * min(ctx$max_abs_NES / 2.5, 1) +
      0.20 * min(rel$n_sig_moa_classes / 10, 1) +
      0.20 * min(rel$n_sig_targets / 10, 1)

    tier <- if (rel$frac_annotated < 0.20 || ctx$n_drugs_sig == 0) {
      "data-limited"
    } else if (score >= 0.70) {
      "strong"
    } else if (score >= 0.45) {
      "moderate"
    } else {
      "weak"
    }

    data.frame(
      contrast = ctx$contrast,
      method = ctx$method,
      score = score,
      tier = tier,
      significant_drugs = ctx$n_drugs_sig,
      max_abs_NES = ctx$max_abs_NES,
      sig_moa = rel$n_sig_moa_classes,
      sig_targets = rel$n_sig_targets,
      annotated_frac = rel$frac_annotated,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out[order(-out$score), , drop = FALSE]
}

drugconnectivity_build_report_tables <- function(pgx,
                                                 method,
                                                 only_annotated = FALSE,
                                                 max_contexts = 8L,
                                                 ntop = 10L) {
  dr <- dc_get_method_data(pgx, method)
  if (is.null(dr)) return(list(text = "No pre-computed drug data available.", data = list(), ranking = data.frame()))

  contrasts <- dc_get_contrasts(dr)
  if (length(contrasts) == 0) return(list(text = "No contrasts available.", data = list(), ranking = data.frame()))

  contexts <- list()
  for (ct in contrasts) {
    ctx <- extract_drugconnectivity_context_data(
      pgx = pgx,
      contrast = ct,
      method = method,
      only_annotated = only_annotated,
      n_top = ntop
    )
    if (!is.null(ctx)) contexts[[ct]] <- ctx
  }

  if (length(contexts) == 0) return(list(text = "No reportable contexts.", data = list(), ranking = data.frame()))

  ranking <- drugconnectivity_rank_contexts(contexts)
  keep <- head(ranking$contrast, max_contexts)
  contexts <- contexts[keep]

  rank_table <- omicsai::omicsai_format_mdtable(
    head(ranking, max_contexts),
    formatters = list(
      score = function(x) omicsai::omicsai_format_num(x, 2),
      max_abs_NES = function(x) omicsai::omicsai_format_num(x, 2),
      annotated_frac = dc_percent
    )
  )

  lines <- c(
    paste0("EXPERIMENT: ", pgx$name %||% pgx$description %||% "omics experiment"),
    paste0("ANALYSIS TYPE: ", dc_analysis_type_info(method)$analysis_type),
    paste0("ANALYSIS TYPE DESCRIPTION: ", dc_analysis_type_info(method)$analysis_type_description),
    "",
    "## Context Ranking",
    rank_table,
    "",
    "## Context Detail"
  )

  for (ct in keep) {
    ctx <- contexts[[ct]]
    rel <- ctx$reliability
    lines <- c(
      lines,
      paste0("### ", ct),
      paste0("- Drugs tested: ", ctx$n_drugs_tested,
             " | significant: ", ctx$n_drugs_sig,
             " | opposing: ", ctx$n_neg,
             " | mimicking: ", ctx$n_pos,
             " | max |NES|: ", omicsai::omicsai_format_num(ctx$max_abs_NES, 2)),
      paste0("- Annotation coverage: ", dc_percent(rel$frac_annotated),
             " | sig MOA classes: ", rel$n_sig_moa_classes,
             " | sig targets: ", rel$n_sig_targets),
      "",
      "Top opposing drugs:",
      dc_compact_drug_table(ctx$top_opposing, n = ntop),
      "",
      "Top mimicking drugs:",
      dc_compact_drug_table(ctx$top_mimicking, n = ntop),
      "",
      "MOA class enrichment:",
      dc_compact_moa_table(ctx$moa_class, n = ntop, label = "MOA Class"),
      "",
      "Target enrichment:",
      dc_compact_moa_table(ctx$moa_target, n = ntop, label = "Target"),
      ""
    )
  }

  list(
    text = paste(lines, collapse = "\n"),
    data = contexts,
    ranking = ranking
  )
}

drugconnectivity_build_methods <- function(pgx, method) {
  template <- omicsai::omicsai_load_template(
    file.path(DRUGCONNECTIVITY_PROMPTS_DIR, "drugconnectivity_report_methods.md")
  )

  params <- list(
    experiment = pgx$name %||% pgx$description %||% "omics experiment",
    analysis_type = dc_analysis_type_info(method)$analysis_type,
    analysis_type_description = dc_analysis_type_info(method)$analysis_type_description,
    date = format(Sys.Date(), "%Y-%m-%d")
  )

  omicsai::omicsai_substitute_template(template, params)
}
