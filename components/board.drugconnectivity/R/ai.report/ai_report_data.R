##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## =============================================================================
## Drug Connectivity AI Report — Data Extraction & Formatting
## -----------------------------------------------------------------------------
## Pure data functions (no Shiny, no LLM). Consumed by ai_text_server.R.
##
## Layered structure:
##   1. Constants and small helpers
##   2. Extraction primitives (call playbase::drugs.* for the four promoted ops)
##   3. Direction subsetting and term selection (board-local; report-specific)
##   4. Context builder (extract_drugconnectivity_context_data)
##   5. Markdown formatters
##   6. Section builders (rank_table, moa_matrix, contrast_block)
##   7. Orchestrators (build_summary_params, build_report_tables, build_methods)
## =============================================================================

DRUGCONNECTIVITY_PROMPTS_DIR <- file.path(OPG, "components/board.drugconnectivity/prompts")

dc_percent <- function(x, digits = 1L) {
  sprintf(paste0("%.", digits, "f%%"), 100 * x)
}

## Internal: split pipe / semicolon / comma separated annotation tokens.
dc_parse_tokens <- function(x) {
  x <- as.character(x %||% "")
  x[is.na(x)] <- ""
  x <- enc2utf8(x)
  lapply(x, function(s) trimws(strsplit(s, split = "[\\|;,]")[[1]]))
}

## -----------------------------------------------------------------------------
## L1000 annotation fallback — reads the bundled CSV when pgx$drugs[[method]]$annot
## is missing. Only used when calling playbase::drugs.enrichmentTable() from the
## board, where FILESX is a board-side constant.
## -----------------------------------------------------------------------------
dc_load_annot_fallback <- function() {
  tryCatch(
    {
      a <- read.csv(
        file.path(FILESX, "cmap/L1000_repurposing_drugs.txt"),
        sep = "\t", comment.char = "#"
      )
      if (!is.null(a) && "pert_iname" %in% colnames(a)) {
        rownames(a) <- a$pert_iname
      }
      a
    },
    error = function(e) NULL
  )
}

## -----------------------------------------------------------------------------
## Direction subset + tier-rank ordering of an MOA / target enrichment table.
## Adds `support` and `support_rank` columns. Used by dc_summary_terms() and
## dc_corroborating_targets().
## -----------------------------------------------------------------------------
dc_enrichment_direction_subset <- function(dt, direction = c("opposing", "mimicking")) {
  direction <- match.arg(direction)
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return(dt[0, , drop = FALSE])
  sel <- if (direction == "opposing") dt$NES < 0 else dt$NES > 0
  out <- dt[sel, , drop = FALSE]
  if (nrow(out) == 0) return(out)
  out$support <- vapply(seq_len(nrow(out)), function(i) {
    playbase::drugs.supportBucket(out$padj[i], out$pval[i])
  }, character(1))
  out$support_rank <- match(out$support, c("significant", "nominal", "unsupported"))
  out[order(out$support_rank, out$padj, out$pval, -abs(out$NES), out$pathway), , drop = FALSE]
}

dc_summary_terms <- function(dt, direction = c("opposing", "mimicking"), n = 3L) {
  dt <- dc_enrichment_direction_subset(dt, direction = direction)
  if (is.null(dt) || nrow(dt) == 0) {
    empty <- if (is.data.frame(dt)) dt[0, , drop = FALSE] else data.frame()
    return(list(
      overall = empty, supported = empty,
      significant = empty, nominal = empty, unsupported = empty
    ))
  }

  significant <- dt[dt$support == "significant", , drop = FALSE]
  nominal     <- dt[dt$support == "nominal", , drop = FALSE]
  unsupported <- dt[dt$support == "unsupported", , drop = FALSE]
  supported <- if (nrow(significant) > 0) {
    significant
  } else if (nrow(nominal) > 0) {
    nominal
  } else {
    dt[0, , drop = FALSE]
  }

  list(
    overall     = head(dt, n),
    supported   = head(supported, n),
    significant = head(significant, n),
    nominal     = head(nominal, n),
    unsupported = head(unsupported, n)
  )
}

dc_corroborating_targets <- function(dsea_table, moa_target, moa_terms,
                                     direction = c("opposing", "mimicking"),
                                     n = 3L) {
  direction <- match.arg(direction)
  if (is.null(dsea_table) || !is.data.frame(dsea_table) || nrow(dsea_table) == 0 ||
      is.null(moa_target) || !is.data.frame(moa_target) || nrow(moa_target) == 0 ||
      is.null(moa_terms) || length(moa_terms) == 0) {
    return(data.frame())
  }

  moa_tokens <- dc_parse_tokens(dsea_table$moa)
  target_tokens <- dc_parse_tokens(dsea_table$target)
  keep <- vapply(moa_tokens, function(tokens) any(tokens %in% moa_terms), logical(1))
  if (!any(keep)) return(moa_target[0, , drop = FALSE])

  supported_targets <- unique(unlist(target_tokens[keep]))
  supported_targets <- setdiff(supported_targets, c("", " ", NA, "NA", "N/A"))
  if (length(supported_targets) == 0) return(moa_target[0, , drop = FALSE])

  dt <- dc_enrichment_direction_subset(moa_target, direction = direction)
  if (is.null(dt) || nrow(dt) == 0) return(dt)
  dt <- dt[dt$pathway %in% supported_targets, , drop = FALSE]
  if (nrow(dt) == 0) return(dt)
  head(dt[order(dt$support_rank, dt$padj, dt$pval, -abs(dt$NES), dt$pathway), , drop = FALSE], n)
}

dc_annotation_confidence <- function(frac_annotated) {
  if (is.na(frac_annotated)) return("unknown")
  if (frac_annotated < 0.20) return("very low")
  if (frac_annotated < 0.30) return("low")
  if (frac_annotated < 0.60) return("moderate")
  "high"
}

dc_build_moa_summary <- function(moa_class, n = 3L) {
  list(
    opposing  = dc_summary_terms(moa_class, direction = "opposing", n = n),
    mimicking = dc_summary_terms(moa_class, direction = "mimicking", n = n)
  )
}

dc_build_target_summary <- function(dsea_table, moa_target, moa_summary, n = 3L) {
  opposing_terms  <- moa_summary$opposing$supported$pathway %||% character(0)
  mimicking_terms <- moa_summary$mimicking$supported$pathway %||% character(0)

  list(
    opposing = list(
      supported = dc_corroborating_targets(
        dsea_table = dsea_table, moa_target = moa_target,
        moa_terms = opposing_terms, direction = "opposing", n = n
      ),
      overall = dc_summary_terms(moa_target, direction = "opposing", n = n)
    ),
    mimicking = list(
      supported = dc_corroborating_targets(
        dsea_table = dsea_table, moa_target = moa_target,
        moa_terms = mimicking_terms, direction = "mimicking", n = n
      ),
      overall = dc_summary_terms(moa_target, direction = "mimicking", n = n)
    )
  )
}

dc_recommended_exemplars <- function(dsea_table, moa_terms,
                                     direction = c("opposing", "mimicking"),
                                     n = 3L) {
  direction <- match.arg(direction)
  if (is.null(dsea_table) || !is.data.frame(dsea_table) || nrow(dsea_table) == 0 ||
      is.null(moa_terms) || length(moa_terms) == 0) {
    return(data.frame())
  }

  dt <- dsea_table[!is.na(dsea_table$NES), , drop = FALSE]
  dt <- if (direction == "opposing") {
    dt[dt$NES < 0, , drop = FALSE]
  } else {
    dt[dt$NES > 0, , drop = FALSE]
  }
  if (nrow(dt) == 0) return(dt)

  moa_tokens <- dc_parse_tokens(dt$moa)
  keep <- vapply(moa_tokens, function(tokens) any(tokens %in% moa_terms), logical(1))
  dt <- dt[keep, , drop = FALSE]
  if (nrow(dt) == 0) return(dt)

  dt <- dt[!(is.na(dt$moa) | dt$moa == ""), , drop = FALSE]
  if (nrow(dt) == 0) return(dt)

  dt$support_rank <- vapply(seq_len(nrow(dt)), function(i) {
    match(
      playbase::drugs.supportBucket(dt$padj[i], dt$pval[i]),
      c("significant", "nominal", "unsupported")
    )
  }, integer(1))

  head(dt[order(dt$support_rank, dt$padj, -abs(dt$NES), dt$drug), , drop = FALSE], n)
}

## -----------------------------------------------------------------------------
## Context builder — extract deterministic data for one contrast/method context.
## -----------------------------------------------------------------------------
extract_drugconnectivity_context_data <- function(pgx, contrast, method,
                                                  dsea_table = NULL,
                                                  moa_class = NULL,
                                                  moa_target = NULL,
                                                  only_annotated = FALSE,
                                                  n_top = 15L) {
  if (is.null(dsea_table)) {
    fallback_annot <- if (is.null(pgx$drugs[[method]]$annot)) dc_load_annot_fallback() else NULL
    dsea_table <- playbase::drugs.enrichmentTable(
      pgx = pgx, method = method, contrast = contrast,
      annot = fallback_annot, only_annotated = only_annotated
    )
  }

  if (is.null(dsea_table) || !is.data.frame(dsea_table) || nrow(dsea_table) == 0) {
    return(NULL)
  }

  if (is.null(moa_class))  moa_class  <- playbase::drugs.moaEnrichment(dsea_table, "moa")
  if (is.null(moa_target)) moa_target <- playbase::drugs.moaEnrichment(dsea_table, "target")

  experiment <- pgx$name %||% pgx$description %||% "omics experiment"
  has_annot <- (dsea_table$moa != "" | dsea_table$target != "")
  has_annot[is.na(has_annot)] <- FALSE

  top_opposing  <- head(dsea_table[order(dsea_table$NES), , drop = FALSE], n_top)
  top_mimicking <- head(dsea_table[order(-dsea_table$NES), , drop = FALSE], n_top)
  top_abs       <- head(dsea_table[order(-abs(dsea_table$NES)), , drop = FALSE], n_top)
  moa_summary   <- dc_build_moa_summary(moa_class, n = 3L)
  target_summary <- dc_build_target_summary(
    dsea_table = dsea_table, moa_target = moa_target,
    moa_summary = moa_summary, n = 3L
  )
  exemplars <- list(
    opposing  = dc_recommended_exemplars(
      dsea_table = dsea_table,
      moa_terms = moa_summary$opposing$supported$pathway %||% character(0),
      direction = "opposing", n = 3L
    ),
    mimicking = dc_recommended_exemplars(
      dsea_table = dsea_table,
      moa_terms = moa_summary$mimicking$supported$pathway %||% character(0),
      direction = "mimicking", n = 3L
    )
  )
  frac_annotated <- mean(has_annot)

  caveats <- c(
    "Connectivity reflects L1000 cell-line perturbation signatures, not direct clinical efficacy.",
    "Dose, cell context, and perturbation time in L1000 may differ from study biology."
  )

  list(
    context_id    = paste(contrast, method, sep = "::"),
    contrast      = contrast,
    method        = method,
    experiment    = experiment,
    n_drugs_tested = nrow(dsea_table),
    n_drugs_sig   = sum(dsea_table$padj < 0.05, na.rm = TRUE),
    n_pos         = sum(dsea_table$NES > 0, na.rm = TRUE),
    n_neg         = sum(dsea_table$NES < 0, na.rm = TRUE),
    max_abs_NES   = max(abs(dsea_table$NES), na.rm = TRUE),
    top_opposing  = top_opposing,
    top_mimicking = top_mimicking,
    top_abs       = top_abs,
    moa_class     = moa_class,
    moa_target    = moa_target,
    moa_summary   = moa_summary,
    target_summary = target_summary,
    exemplars     = exemplars,
    reliability = list(
      has_annotations = any(has_annot),
      frac_annotated  = frac_annotated,
      annotation_confidence = dc_annotation_confidence(frac_annotated),
      n_sig_moa_classes = if (is.data.frame(moa_class)) sum(moa_class$padj < 0.05, na.rm = TRUE) else 0L,
      n_sig_targets     = if (is.data.frame(moa_target)) sum(moa_target$padj < 0.05, na.rm = TRUE) else 0L
    ),
    caveats = caveats
  )
}

## -----------------------------------------------------------------------------
## Markdown formatters — render context fields into compact tables and prose.
## -----------------------------------------------------------------------------
dc_compact_drug_table <- function(dt, n = 10L) {
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return("No data available.")
  dt <- head(dt, n)
  fmt <- data.frame(
    Drug = dt$drug,
    NES = dt$NES,
    `q-value` = dt$padj,
    MOA = ifelse(is.na(dt$moa) | dt$moa == "", "-", dt$moa),
    Target = ifelse(is.na(dt$target) | dt$target == "", "-", dt$target),
    check.names = FALSE, stringsAsFactors = FALSE
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
    check.names = FALSE, stringsAsFactors = FALSE
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
    paste0("- Supported opposing MOA terms: ",
           dc_format_summary_entries(moa_summary$opposing$supported, n = 3L)),
    paste0("- Supported mimicking MOA terms: ",
           dc_format_summary_entries(moa_summary$mimicking$supported, n = 3L)),
    paste0("- Corroborating opposing targets: ",
           dc_format_summary_entries(target_summary$opposing$supported, n = 3L)),
    paste0("- Corroborating mimicking targets: ",
           dc_format_summary_entries(target_summary$mimicking$supported, n = 3L)),
    paste0("- Preferred opposing exemplars matching supported MOA terms: ",
           dc_format_exemplars(ctx$exemplars$opposing, n = 3L)),
    paste0("- Preferred mimicking exemplars matching supported MOA terms: ",
           dc_format_exemplars(ctx$exemplars$mimicking, n = 3L)),
    paste0("- Annotation confidence: ", rel$annotation_confidence,
           " (", dc_percent(rel$frac_annotated), " annotated)")
  ), collapse = "\n")
}

## -----------------------------------------------------------------------------
## Render a single per-contrast block from the shared inner template
## (drugconnectivity_contrast_data.md). Used by both Summary mode (one block)
## and Report mode (one block per contrast inside the multi-contrast scaffold).
## -----------------------------------------------------------------------------
dc_render_contrast_block <- function(ctx, tier_label, ntop = 10L) {
  if (is.null(ctx)) return("")
  rel <- ctx$reliability

  headline_metrics <- paste0(
    "Drugs tested: ", ctx$n_drugs_tested,
    " | significant (q<0.05): ", ctx$n_drugs_sig,
    " | opposing (neg NES): ", ctx$n_neg,
    " | mimicking (pos NES): ", ctx$n_pos,
    " | max |NES|: ", omicsai::omicsai_format_num(ctx$max_abs_NES, 2)
  )

  annotation_coverage <- paste0(
    dc_percent(rel$frac_annotated),
    " | sig MOA classes: ", rel$n_sig_moa_classes,
    " | sig targets: ", rel$n_sig_targets
  )

  template <- omicsai::omicsai_load_template(
    file.path(DRUGCONNECTIVITY_PROMPTS_DIR, "drugconnectivity_contrast_data.md")
  )

  omicsai::omicsai_substitute_template(template, list(
    contrast            = ctx$contrast,
    tier                = tier_label,
    headline_metrics    = headline_metrics,
    annotation_coverage = annotation_coverage,
    evidence_summary    = dc_evidence_summary_block(ctx),
    top_opposing_table  = dc_compact_drug_table(ctx$top_opposing, n = ntop),
    top_mimicking_table = dc_compact_drug_table(ctx$top_mimicking, n = ntop),
    moa_class_table     = dc_compact_moa_table(ctx$moa_class, n = ntop, label = "MOA Class"),
    moa_target_table    = dc_compact_moa_table(ctx$moa_target, n = ntop, label = "Target")
  ))
}

## -----------------------------------------------------------------------------
## Cross-contrast MOA convergence matrix.
## Rows = top MOA classes (ranked by max |NES| across contrasts).
## Columns = contrasts.  Cells = NES with significance marker (** q<0.05, * p<0.05).
## -----------------------------------------------------------------------------
dc_moa_convergence_matrix <- function(contexts, n_moa = 12L) {
  if (length(contexts) == 0) return(NULL)

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
    if (is.null(entry) || is.null(entry$NES) || is.na(entry$NES)) return("—")
    sig <- if (!is.na(entry$padj) && entry$padj < 0.05) "**"
           else if (!is.na(entry$pval) && entry$pval < 0.05) "*"
           else ""
    paste0(omicsai::omicsai_format_num(entry$NES, 2), sig)
  }

  col_names  <- c("MOA Class", contrasts)
  header_row <- paste0("| ", paste(col_names, collapse = " | "), " |")
  sep_row    <- paste0("| ", paste(rep("---", length(col_names)), collapse = " | "), " |")
  data_rows  <- vapply(top_moa, function(moa) {
    cells <- vapply(contrasts, function(ct) fmt_cell(moa_entries[[moa]][[ct]]), character(1))
    paste0("| ", paste(c(moa, cells), collapse = " | "), " |")
  }, character(1))

  paste(c(header_row, sep_row, data_rows), collapse = "\n")
}

## -----------------------------------------------------------------------------
## Orchestrators — assemble template-ready params for Summary and Report modes.
## -----------------------------------------------------------------------------
drugconnectivity_build_summary_params <- function(pgx, contrast, method,
                                                  only_annotated = FALSE,
                                                  ntop = 12L) {
  ctx <- extract_drugconnectivity_context_data(
    pgx = pgx, contrast = contrast, method = method,
    only_annotated = only_annotated, n_top = ntop
  )
  if (is.null(ctx)) return(NULL)

  ranking <- drugconnectivity_rank_contexts(list(default = ctx))
  tier_label <- if (nrow(ranking) > 0) ranking$tier[1] else "unknown"

  at <- playbase::drugs.analysisInfo(method)

  list(
    experiment                = ctx$experiment,
    analysis_type             = at$analysis_type,
    analysis_type_description = at$analysis_type_description,
    contrast_detail           = dc_render_contrast_block(ctx, tier_label, ntop = ntop)
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

drugconnectivity_build_report_tables <- function(pgx, method,
                                                 only_annotated = FALSE,
                                                 max_contexts = 8L,
                                                 ntop = 10L) {
  dr <- pgx$drugs[[method]]
  if (is.null(dr)) {
    return(list(text = "No pre-computed drug data available.", data = list(), ranking = data.frame()))
  }

  contrasts <- if (!is.null(dr$X)) sort(setdiff(colnames(dr$X), grep("^IA:", colnames(dr$X), value = TRUE))) else character(0)
  if (length(contrasts) == 0) {
    return(list(text = "No contrasts available.", data = list(), ranking = data.frame()))
  }

  contexts <- list()
  for (ct in contrasts) {
    ctx <- extract_drugconnectivity_context_data(
      pgx = pgx, contrast = ct, method = method,
      only_annotated = only_annotated, n_top = ntop
    )
    if (!is.null(ctx)) contexts[[ct]] <- ctx
  }

  if (length(contexts) == 0) {
    return(list(text = "No reportable contexts.", data = list(), ranking = data.frame()))
  }

  ranking <- drugconnectivity_rank_contexts(contexts)
  keep <- head(ranking$contrast, max_contexts)
  contexts <- contexts[keep]

  rank_table <- paste(omicsai::omicsai_format_mdtable(
    head(ranking, max_contexts),
    formatters = list(
      score        = function(x) omicsai::omicsai_format_num(x, 2),
      max_abs_NES  = function(x) omicsai::omicsai_format_num(x, 2),
      annotated_frac = dc_percent
    )
  ), collapse = "\n")

  moa_matrix <- dc_moa_convergence_matrix(contexts, n_moa = 12L) %||% "(no MOA convergence data)"

  at <- playbase::drugs.analysisInfo(method)

  contrast_blocks <- vapply(keep, function(ct) {
    tier_row <- ranking[ranking$contrast == ct, , drop = FALSE]
    tier_label <- if (nrow(tier_row) > 0) tier_row$tier[1] else "unknown"
    dc_render_contrast_block(contexts[[ct]], tier_label, ntop = ntop)
  }, character(1))

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

  list(text = text, data = contexts, ranking = ranking)
}

drugconnectivity_build_methods <- function(pgx, method) {
  template <- omicsai::omicsai_load_template(
    file.path(DRUGCONNECTIVITY_PROMPTS_DIR, "drugconnectivity_methods.md")
  )

  at <- playbase::drugs.analysisInfo(method)
  params <- list(
    experiment = pgx$name %||% pgx$description %||% "omics experiment",
    analysis_type = at$analysis_type,
    analysis_type_description = at$analysis_type_description,
    date = format(Sys.Date(), "%Y-%m-%d")
  )

  omicsai::omicsai_substitute_template(template, params)
}
