##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

DRUGCONNECTIVITY_PROMPTS_DIR <- file.path(OPG, "components/board.drugconnectivity/prompts")

dc_percent <- function(x, digits = 1L) {
  sprintf(paste0("%.", digits, "f%%"), 100 * x)
}

dc_analysis_type_info <- function(analysis_type) {
  # Normalize to lowercase for case-insensitive matching
  key <- tolower(analysis_type %||% "")

  desc <- if (key %in% c("l1000/activity", "l1000_activitys_n20d1011")) {
    paste(
      "L1000 Transcriptional Activity Scores (TAS) connectivity.",
      "Each drug signature is weighted by strength and reproducibility across cell lines.",
      "Negative NES = drug opposes the experimental state (reversal candidate);",
      "positive NES = drug mimics it (shared transcriptional programme)."
    )
  } else if (key == "l1000/gene") {
    paste(
      "L1000 gene-level log-fold-change connectivity (978 landmark genes, per-cell-line).",
      "More granular than TAS but noisier; best used to corroborate TAS findings.",
      "Negative NES = drug opposes the experimental state;",
      "positive NES = drug mimics it."
    )
  } else if (key == "ctrpv2/sensitivity") {
    paste(
      "CTRPv2 pharmacological sensitivity connectivity (AUC dose-response, ~481 compounds, ~860 cancer cell lines).",
      "NES reflects co-variation of drug sensitivity with the experimental signature.",
      "Positive NES = the experimental state predicts drug vulnerability (sensitivity);",
      "negative NES = resistance or insensitivity."
    )
  } else if (key == "gdsc/sensitivity") {
    paste(
      "GDSC pharmacogenomic sensitivity connectivity (IC50-based, ~367 clinical anti-cancer drugs, ~987 cell lines).",
      "Strongest clinical relevance due to approved and late-stage compounds.",
      "Positive NES = the experimental state predicts drug sensitivity;",
      "negative NES = predicted resistance."
    )
  } else {
    paste(
      "Drug connectivity analysis based on the selected pre-computed signature resource.",
      "Interpret results according to the selected analysis type metadata."
    )
  }

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
  paste(omicsai::omicsai_format_mdtable(
    fmt,
    formatters = list(
      NES = function(x) omicsai::omicsai_format_num(x, 3),
      `q-value` = omicsai::omicsai_format_pvalue
    )
  ), collapse = "\n")
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
  paste(omicsai::omicsai_format_mdtable(
    fmt,
    formatters = list(
      NES = function(x) omicsai::omicsai_format_num(x, 3),
      `q-value` = omicsai::omicsai_format_pvalue
    )
  ), collapse = "\n")
}

dc_summary_label <- function(dt) {
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return("no supported MOA term")
  label <- dt$pathway[1]
  support <- dt$support[1] %||% "unsupported"
  if (identical(support, "significant")) return(label)
  paste0(label, " (", support, ")")
}

dc_format_summary_entries <- function(dt, n = 3L, include_support = TRUE) {
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return("none")
  dt <- head(dt, n)
  entries <- vapply(seq_len(nrow(dt)), function(i) {
    parts <- c(
      dt$pathway[i],
      paste0("NES=", omicsai::omicsai_format_num(dt$NES[i], 2))
    )
    if (!is.na(dt$padj[i])) {
      parts <- c(parts, paste0("q=", omicsai::omicsai_format_pvalue(dt$padj[i])))
    }
    if (isTRUE(include_support) && "support" %in% colnames(dt)) {
      parts <- c(parts, dt$support[i])
    }
    paste0(parts[1], " (", paste(parts[-1], collapse = ", "), ")")
  }, character(1))
  paste(entries, collapse = "; ")
}

dc_format_exemplars <- function(dt, n = 3L) {
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return("none")
  dt <- head(dt, n)
  entries <- vapply(seq_len(nrow(dt)), function(i) {
    parts <- c(
      dt$drug[i],
      dt$moa[i],
      paste0("NES=", omicsai::omicsai_format_num(dt$NES[i], 2))
    )
    if (!is.na(dt$padj[i])) {
      parts <- c(parts, paste0("q=", omicsai::omicsai_format_pvalue(dt$padj[i])))
    }
    paste0(parts[1], " [", paste(parts[-1], collapse = ", "), "]")
  }, character(1))
  paste(entries, collapse = "; ")
}

dc_evidence_summary_block <- function(ctx) {
  rel <- ctx$reliability
  moa_summary <- ctx$moa_summary
  target_summary <- ctx$target_summary

  paste(c(
    "Interpretation evidence summary (primary evidence for report writing):",
    paste0(
      "- Supported opposing MOA terms: ",
      dc_format_summary_entries(moa_summary$opposing$supported, n = 3L)
    ),
    paste0(
      "- Supported mimicking MOA terms: ",
      dc_format_summary_entries(moa_summary$mimicking$supported, n = 3L)
    ),
    paste0(
      "- Corroborating opposing targets: ",
      dc_format_summary_entries(target_summary$opposing$supported, n = 3L)
    ),
    paste0(
      "- Corroborating mimicking targets: ",
      dc_format_summary_entries(target_summary$mimicking$supported, n = 3L)
    ),
    paste0(
      "- Preferred opposing exemplars matching supported MOA terms: ",
      dc_format_exemplars(ctx$exemplars$opposing, n = 3L)
    ),
    paste0(
      "- Preferred mimicking exemplars matching supported MOA terms: ",
      dc_format_exemplars(ctx$exemplars$mimicking, n = 3L)
    ),
    paste0(
      "- Annotation confidence: ", rel$annotation_confidence,
      " (", dc_percent(rel$frac_annotated), " annotated)"
    )
  ), collapse = "\n")
}

# One-liner summary for a single contrast — used as section sub-header in data block.
# Identifies the dominant opposing and mimicking MOA theme from moa_class enrichment.
dc_contrast_oneliner <- function(ctx, tier) {
  n_sig   <- ctx$n_drugs_sig %||% 0L
  max_nes <- ctx$max_abs_NES %||% 0

  opp_theme <- dc_summary_label(ctx$moa_summary$opposing$supported)
  mim_theme <- dc_summary_label(ctx$moa_summary$mimicking$supported)

  paste0(
    tier, " | ", n_sig, " sig drugs | max|NES|=",
    omicsai::omicsai_format_num(max_nes, 2),
    " | top opposing: ", opp_theme,
    " | top mimicking: ", mim_theme
  )
}

# Cross-contrast MOA convergence matrix.
# Rows = top MOA classes (ranked by max |NES| across contrasts).
# Columns = contrasts.  Cells = NES with significance marker (** q<0.05, * p<0.05).
# Generalises to any number of contrasts (works for 1 or many).
dc_moa_convergence_matrix <- function(contexts, n_moa = 12L) {
  if (length(contexts) == 0) return(NULL)

  # Collect MOA enrichment: moa_name -> list(contrast_name -> list(NES, padj, pval))
  moa_entries <- list()
  for (ct in names(contexts)) {
    mc <- contexts[[ct]]$moa_class
    if (is.null(mc) || nrow(mc) == 0) next
    for (i in seq_len(nrow(mc))) {
      moa <- mc$pathway[i]
      if (!moa %in% names(moa_entries)) moa_entries[[moa]] <- list()
      moa_entries[[moa]][[ct]] <- list(
        NES  = mc$NES[i],
        padj = mc$padj[i],
        pval = mc$pval[i]
      )
    }
  }
  if (length(moa_entries) == 0) return(NULL)

  # Rank MOA classes by recurrence of supported evidence, then max |NES|.
  scores <- vapply(moa_entries, function(entries) {
    sig_n <- sum(vapply(entries, function(e) !is.na(e$padj) && e$padj < 0.05, logical(1)))
    nom_n <- sum(vapply(entries, function(e) {
      !is.na(e$padj) && e$padj >= 0.05 && !is.na(e$pval) && e$pval < 0.05
    }, logical(1)))
    max_abs <- max(abs(vapply(entries, function(e) e$NES %||% 0, numeric(1))), na.rm = TRUE)
    100 * sig_n + 10 * nom_n + max_abs
  }, numeric(1))
  top_moa <- names(sort(scores, decreasing = TRUE))[seq_len(min(n_moa, length(scores)))]
  contrasts <- names(contexts)

  fmt_cell <- function(entry) {
    if (is.null(entry) || is.null(entry$NES) || is.na(entry$NES)) return("\u2014")
    sig <- if (!is.na(entry$padj) && entry$padj < 0.05) "**"
           else if (!is.na(entry$pval) && entry$pval < 0.05) "*"
           else ""
    paste0(omicsai::omicsai_format_num(entry$NES, 2), sig)
  }

  # Build markdown table as raw strings (values are pre-formatted)
  col_names  <- c("MOA Class", contrasts)
  header_row <- paste0("| ", paste(col_names, collapse = " | "), " |")
  sep_row    <- paste0("| ", paste(rep("---", length(col_names)), collapse = " | "), " |")
  data_rows  <- vapply(top_moa, function(moa) {
    cells <- vapply(contrasts, function(ct) fmt_cell(moa_entries[[moa]][[ct]]), character(1))
    paste0("| ", paste(c(moa, cells), collapse = " | "), " |")
  }, character(1))

  paste(c(header_row, sep_row, data_rows), collapse = "\n")
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

  # --- Computed parts (R handles numbers and tables; .md handles labels) ---

  rank_table <- paste(omicsai::omicsai_format_mdtable(
    head(ranking, max_contexts),
    formatters = list(
      score        = function(x) omicsai::omicsai_format_num(x, 2),
      max_abs_NES  = function(x) omicsai::omicsai_format_num(x, 2),
      annotated_frac = dc_percent
    )
  ), collapse = "\n")

  moa_matrix <- dc_moa_convergence_matrix(contexts, n_moa = 12L) %||% "(no MOA convergence data)"

  at <- dc_analysis_type_info(method)

  # --- Per-contrast blocks (repeating unit — built in R) ---
  contrast_blocks <- vapply(keep, function(ct) {
    ctx      <- contexts[[ct]]
    rel      <- ctx$reliability
    tier_row <- ranking[ranking$contrast == ct, , drop = FALSE]
    tier_label <- if (nrow(tier_row) > 0) tier_row$tier[1] else "unknown"

    paste(c(
      paste0("### ", ct),
      paste0("(", dc_contrast_oneliner(ctx, tier_label), ")"),
      paste0("- Drugs tested: ", ctx$n_drugs_tested,
             " | significant (q<0.05): ", ctx$n_drugs_sig,
             " | opposing (neg NES): ", ctx$n_neg,
             " | mimicking (pos NES): ", ctx$n_pos,
             " | max |NES|: ", omicsai::omicsai_format_num(ctx$max_abs_NES, 2)),
      paste0("- Annotation coverage: ", dc_percent(rel$frac_annotated),
             " | sig MOA classes: ", rel$n_sig_moa_classes,
             " | sig targets: ", rel$n_sig_targets),
      "",
      dc_evidence_summary_block(ctx),
      "",
      "Top opposing drugs (negative NES \u2014 candidate signature reversal):",
      dc_compact_drug_table(ctx$top_opposing, n = ntop),
      "",
      "Top mimicking drugs (positive NES \u2014 mechanistic similarity to experimental state):",
      dc_compact_drug_table(ctx$top_mimicking, n = ntop),
      "",
      "MOA class enrichment (which drug mechanisms converge on this signature):",
      dc_compact_moa_table(ctx$moa_class, n = ntop, label = "MOA Class"),
      "",
      "Target enrichment (which molecular targets drive the connectivity):",
      dc_compact_moa_table(ctx$moa_target, n = ntop, label = "Target"),
      ""
    ), collapse = "\n")
  }, character(1))

  # --- Full data block — rendered from single report_data template ---
  header_tmpl <- omicsai::omicsai_load_template(
    file.path(DRUGCONNECTIVITY_PROMPTS_DIR, "drugconnectivity_report_data.md")
  )

  text <- omicsai::omicsai_substitute_template(header_tmpl, list(
    experiment               = pgx$name %||% pgx$description %||% "omics experiment",
    analysis_type            = at$analysis_type,
    analysis_type_description = at$analysis_type_description,
    rank_table               = rank_table,
    moa_matrix               = moa_matrix,
    contrast_detail          = paste(contrast_blocks, collapse = "\n")
  ))

  list(
    text    = text,
    data    = contexts,
    ranking = ranking
  )
}

drugconnectivity_build_methods <- function(pgx, method) {
  template <- omicsai::omicsai_load_template(
    file.path(DRUGCONNECTIVITY_PROMPTS_DIR, "drugconnectivity_methods.md")
  )

  params <- list(
    experiment = pgx$name %||% pgx$description %||% "omics experiment",
    analysis_type = dc_analysis_type_info(method)$analysis_type,
    analysis_type_description = dc_analysis_type_info(method)$analysis_type_description,
    date = format(Sys.Date(), "%Y-%m-%d")
  )

  omicsai::omicsai_substitute_template(template, params)
}
