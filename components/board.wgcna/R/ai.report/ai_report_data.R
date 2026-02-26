## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.

# =============================================================================
# WGCNA AI Report — Data Extraction and Classification
# =============================================================================
# Deterministic functions that build structured data for the AI report prompt.
# No Shiny dependencies. No LLM calls.
#
# Helpers (in ai_report_data_extract.R):
#   extract_module_data, resolve_symbols, resolve_functions,
#   select_top_enrichment, classify_artifact
#
# Public functions:
#   wgcna_build_report_tables(wgcna, pgx, ...)  — structured data tables
#   wgcna_build_summary_params(wgcna, module, pgx) — single-module summary params
#   wgcna_rank_modules(wgcna)                    — signal tier classification
#   wgcna_build_methods(wgcna, pgx)              — deterministic methods section

BOARD_PROMPTS_DIR <- file.path(OPG, "components/board.wgcna/prompts")

#' Build structured report tables from WGCNA results
#'
#' Extracts module-level data (enrichment, hub genes, eigengene profiles,
#' trait correlations) and formats them as a text block for the LLM prompt
#' plus a structured list for deterministic use.
#'
#' @param wgcna WGCNA results object
#' @param pgx PGX object
#' @param n_modules Integer; max non-grey modules to detail (default 8)
#' @param ntop_enrichment Integer; max enrichment terms per module (default 5)
#' @param ntop_genes Integer; max hub genes per module (default 5)
#' @param include_module_cors Logical; include module-module correlations
#' @param include_overlap Logical; include enrichment overlap details
#' @param include_families Logical; include gene family enrichment
#' @param include_contrasts Logical; include experimental contrasts
#'
#' @return List with \code{text} (character) and \code{data} (structured list)
wgcna_build_report_tables <- function(wgcna, pgx,
                                      n_modules = 8L,
                                      ntop_enrichment = 5L,
                                      ntop_genes = 5L,
                                      include_module_cors = TRUE,
                                      include_overlap = TRUE,
                                      include_families = TRUE,
                                      include_contrasts = TRUE) {
  annot <- wgcna$annot

  # --- Experiment metadata ---
  experiment <- wgcna$experiment %||% pgx$description %||% "omics experiment"
  organism <- pgx$organism %||% "unknown"
  n_samples <- tryCatch(nrow(pgx$samples), error = function(e) NA_integer_)
  sample_groups <- if (!is.null(pgx$samples$group)) {
    levels_or_unique <- if (is.factor(pgx$samples$group)) {
      levels(pgx$samples$group)
    } else {
      unique(pgx$samples$group)
    }
    paste(levels_or_unique, collapse = ", ")
  } else {
    "unknown"
  }
  n_features <- tryCatch(nrow(pgx$X), error = function(e) NA_integer_)
  n_wgcna_features <- length(unlist(wgcna$me.genes))
  if (is.null(n_wgcna_features) || n_wgcna_features == 0) {
    n_wgcna_features <- tryCatch(ncol(wgcna$datExpr), error = function(e) NA_integer_)
  }

  # WGCNA parameters
  network_type <- wgcna$net$networkType %||% wgcna$networkType %||% "signed"
  power <- wgcna$power %||% wgcna$net$power %||% "NA"
  min_mod_size <- wgcna$minModSize %||% "20"
  merge_cut_height <- wgcna$mergeCutHeight %||% "0.15"

  # --- Module ordering ---
  all_modules <- names(wgcna$me.genes)
  non_grey <- setdiff(all_modules, "MEgrey")
  non_grey <- non_grey[order(lengths(wgcna$me.genes[non_grey]), decreasing = TRUE)]
  non_grey <- head(non_grey, n_modules)

  # --- Pre-compute shared data for loop efficiency ---
  M <- tryCatch(playbase:::wgcna.get_modTraits(wgcna), error = function(e) NULL)
  top <- tryCatch(
    playbase::wgcna.getTopGenesAndSets(wgcna, annot = annot, ntop = 40,
                                       level = "gene", rename = "gene_title"),
    error = function(e) list(pheno = list(), genes = list(), sets = list())
  )

  # --- Build per-module data ---
  module_data <- list()
  overview_rows <- list()

  for (mod in non_grey) {
    md <- extract_module_data(wgcna, mod, pgx, top = top, M = M)

    # Tiered enrichment selection for report
    top_enrich <- NULL
    if (is.data.frame(md$gse) && md$n_sig > 0) {
      sig_gse <- md$gse[!is.na(md$gse$q.value) & md$gse$q.value < 0.05, , drop = FALSE]
      sig_gse <- sig_gse[order(sig_gse$q.value), , drop = FALSE]
      top_enrich <- select_top_enrichment(sig_gse, ntop_enrichment)
    }

    # Hub genes: use first available trait for single-trait report view
    hub_genes_df <- NULL
    trait_for_genes <- if (length(md$traits) > 0) md$traits[1] else md$top_trait

    gs <- if (nzchar(trait_for_genes) && trait_for_genes %in% names(md$trait_stats)) {
      md$trait_stats[[trait_for_genes]]
    } else {
      NULL
    }

    if (is.data.frame(gs) && nrow(gs) > 0) {
      if ("moduleMembership" %in% colnames(gs)) {
        gs <- gs[order(-abs(gs$moduleMembership)), , drop = FALSE]
      }
      top_gs <- head(gs, ntop_genes)

      symbols <- resolve_symbols(top_gs$feature, annot)
      funcs <- resolve_functions(top_gs$feature, annot)

      mm_vals <- if ("moduleMembership" %in% colnames(top_gs)) {
        top_gs$moduleMembership
      } else {
        rep(NA_real_, nrow(top_gs))
      }
      logfc_vals <- if ("foldChange" %in% colnames(top_gs)) {
        top_gs$foldChange
      } else {
        rep(NA_real_, nrow(top_gs))
      }

      hub_genes_df <- data.frame(
        symbol = symbols,
        MM = mm_vals,
        logFC = logfc_vals,
        func = funcs,
        stringsAsFactors = FALSE
      )
    }

    module_data[[mod]] <- list(
      size = md$size,
      eigengene_profile = md$eigengene_profile,
      peak_condition = md$peak_condition,
      top_trait = md$top_trait,
      top_r = md$top_r,
      n_sig = md$n_sig,
      n_total = md$n_total,
      top_enrichment = top_enrich,
      hub_genes = hub_genes_df,
      trait_for_genes = trait_for_genes
    )

    overview_rows[[mod]] <- list(
      module = mod,
      size = md$size,
      peak_condition = md$peak_condition,
      top_trait = md$top_trait,
      top_r = md$top_r,
      n_sig = md$n_sig,
      n_total = md$n_total
    )
  }

  # --- Contrasts (experimental design) ---
  contrasts_text <- NULL
  if (include_contrasts && !is.null(pgx$model.parameters$contr.matrix)) {
    cm <- pgx$model.parameters$contr.matrix
    contrast_names <- colnames(cm)
    group_names <- rownames(cm)
    ct_lines <- character(0)
    for (cn in contrast_names) {
      pos <- group_names[cm[, cn] > 0]
      neg <- group_names[cm[, cn] < 0]
      ct_lines <- c(ct_lines, sprintf("- %s: %s vs %s",
        cn, paste(pos, collapse = "+"), paste(neg, collapse = "+")))
    }
    contrasts_text <- ct_lines
  }

  # --- Module-module eigengene correlations ---
  module_cors_text <- NULL
  if (include_module_cors) {
    me_data <- wgcna$datME %||% wgcna$net$MEs
    if (!is.null(me_data)) {
      me_cols <- intersect(non_grey, colnames(me_data))
      if (length(me_cols) >= 2) {
        me_cor <- cor(me_data[, me_cols])
        pairs <- character(0)
        for (i in 1:(length(me_cols) - 1)) {
          for (j in (i + 1):length(me_cols)) {
            r <- me_cor[i, j]
            if (abs(r) >= 0.7) {
              pairs <- c(pairs, sprintf("%s \u2194 %s: r = %+.2f",
                me_cols[i], me_cols[j], r))
            }
          }
        }
        module_cors_text <- pairs
      }
    }
  }

  # --- Enrichment overlap details ---
  if (include_overlap) {
    for (mod in non_grey) {
      mdat <- module_data[[mod]]
      if (!is.null(mdat$top_enrichment) && nrow(mdat$top_enrichment) > 0) {
        te <- mdat$top_enrichment
        if ("overlap" %in% colnames(te)) {
          module_data[[mod]]$enrichment_overlap <- te$overlap
        }
        if ("genes" %in% colnames(te)) {
          gene_lists <- vapply(te$genes, function(g) {
            genes <- strsplit(as.character(g), "\\|")[[1]]
            if (length(genes) > 10) {
              paste0(paste(head(genes, 10), collapse = "|"),
                     sprintf(" (+%d more)", length(genes) - 10))
            } else {
              paste(genes, collapse = "|")
            }
          }, character(1))
          module_data[[mod]]$enrichment_genes <- gene_lists
        }
      }
    }
  }

  # --- Gene family enrichment per module ---
  families_text <- NULL
  if (include_families && !is.null(pgx$families)) {
    fam_names <- names(pgx$families)
    fam_names <- fam_names[fam_names != "<all>"]
    fam_sizes <- lengths(pgx$families[fam_names])
    fam_names <- fam_names[fam_sizes >= 5 & fam_sizes <= 500]

    families_text <- list()
    for (mod in non_grey) {
      mod_genes <- wgcna$me.genes[[mod]]
      if (is.null(mod_genes) || length(mod_genes) == 0) next

      overlaps <- vapply(fam_names, function(fn) {
        length(intersect(mod_genes, pgx$families[[fn]]))
      }, integer(1))

      sig_fam <- overlaps[overlaps >= 3]
      if (length(sig_fam) == 0) next
      sig_fam <- sort(sig_fam, decreasing = TRUE)
      sig_fam <- head(sig_fam, 5)

      fam_strs <- vapply(names(sig_fam), function(fn) {
        sprintf("%s (%d of %d)", fn, sig_fam[fn], length(pgx$families[[fn]]))
      }, character(1))

      families_text[[mod]] <- fam_strs
    }
  }

  # --- Format text output ---
  lines <- character(0)

  # Header
  lines <- c(lines,
    paste0("EXPERIMENT: ", experiment),
    paste0("ORGANISM: ", organism),
    paste0("SAMPLES: ", n_samples, " (", sample_groups, ")"),
    paste0("FEATURES: ", n_features, " (", n_wgcna_features, " used for WGCNA)"),
    paste0("WGCNA PARAMETERS: ", network_type,
           ", power=", power,
           ", minModSize=", min_mod_size,
           ", mergeCutHeight=", merge_cut_height),
    ""
  )

  # Overview table sorted by sig enrichments descending
  overview_order <- names(sort(vapply(overview_rows, function(x) x$n_sig, integer(1)),
                               decreasing = TRUE))

  grey_size <- length(wgcna$me.genes[["MEgrey"]])
  ov_rows <- lapply(overview_order, function(m) {
    ov <- overview_rows[[m]]
    data.frame(
      Module = m,
      Genes = as.character(ov$size),
      Peak_Condition = ov$peak_condition,
      Top_Trait = ov$top_trait,
      Trait_r = if (is.na(ov$top_r)) "NA" else sprintf("%+.2f", ov$top_r),
      Sig_Enrichments = as.character(ov$n_sig),
      Total = as.character(ov$n_total),
      stringsAsFactors = FALSE
    )
  })
  ov_rows[[length(ov_rows) + 1]] <- data.frame(
    Module = "grey", Genes = as.character(grey_size),
    Peak_Condition = "\u2014", Top_Trait = "\u2014", Trait_r = "\u2014",
    Sig_Enrichments = "\u2014", Total = "\u2014",
    stringsAsFactors = FALSE
  )
  overview_df <- do.call(rbind, ov_rows)
  lines <- c(lines, "MODULE OVERVIEW:")
  lines <- c(lines, omicsai::omicsai_format_mdtable(overview_df))
  lines <- c(lines, "")

  # Per-module detail
  lines <- c(lines, "PER-MODULE DETAIL:")

  for (mod in overview_order) {
    mdat <- module_data[[mod]]
    lines <- c(lines, sprintf("### %s (%d genes)", mod, mdat$size))

    # Eigengene profile
    if (!is.null(mdat$eigengene_profile)) {
      profile_str <- paste(
        sprintf("%s=%s", names(mdat$eigengene_profile),
                omicsai::omicsai_format_num(mdat$eigengene_profile, 2)),
        collapse = ", "
      )
      lines <- c(lines, paste0("Eigengene profile: ", profile_str))
    }

    # Top trait
    if (nzchar(mdat$top_trait) && !is.na(mdat$top_r)) {
      lines <- c(lines, sprintf("Top trait: %s (r=%+.2f)", mdat$top_trait, mdat$top_r))
    }
    lines <- c(lines, "")

    # Enrichment
    if (!is.null(mdat$top_enrichment) && nrow(mdat$top_enrichment) > 0) {
      te <- mdat$top_enrichment
      enrich_df <- data.frame(
        `#` = seq_len(nrow(te)),
        Term = te$geneset,
        `q-value` = te$q.value,
        Source = if ("source" %in% colnames(te)) as.character(te$source) else "",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      lines <- c(lines,
        sprintf("Top enrichment (%d significant of %d):", mdat$n_sig, mdat$n_total),
        omicsai::omicsai_format_mdtable(enrich_df, formatters = list(
          `q-value` = omicsai::omicsai_format_pvalue
        ))
      )
    } else {
      lines <- c(lines, "No significant enrichment (all q > 0.05)")
    }
    lines <- c(lines, "")

    # Hub genes
    if (!is.null(mdat$hub_genes) && nrow(mdat$hub_genes) > 0) {
      hg <- mdat$hub_genes
      logfc_label <- if (nzchar(mdat$trait_for_genes)) {
        paste0("logFC (vs ", mdat$trait_for_genes, ")")
      } else {
        "logFC"
      }
      lines <- c(lines,
        sprintf("Hub genes (top %d by MM):", nrow(hg)),
        omicsai::omicsai_format_mdtable(
          hg,
          col_labels = c(symbol = "Gene", logFC = logfc_label, func = "Known function"),
          formatters = list(
            MM = function(x) omicsai::omicsai_format_num(x, 2),
            logFC = function(x) omicsai::omicsai_format_num(x, 1)
          )
        )
      )
    }

    # Enrichment overlap (within per-module block)
    if (include_overlap && !is.null(module_data[[mod]]$enrichment_overlap)) {
      ov <- module_data[[mod]]$enrichment_overlap
      gl <- module_data[[mod]]$enrichment_genes
      if (!is.null(ov) && length(ov) > 0) {
        lines <- c(lines, "Enrichment overlap:")
        te <- module_data[[mod]]$top_enrichment
        for (i in seq_along(ov)) {
          line <- sprintf("  %s: %s", te$geneset[i], ov[i])
          if (!is.null(gl) && length(gl) >= i) {
            line <- paste0(line, " \u2014 genes: ", gl[i])
          }
          lines <- c(lines, line)
        }
      }
    }

    # Gene families (within per-module block)
    if (include_families && !is.null(families_text[[mod]])) {
      lines <- c(lines, sprintf("Gene family enrichment: %s",
                                paste(families_text[[mod]], collapse = "; ")))
    }

    lines <- c(lines, "")
  }

  # --- Optional sections (NULL-safe via collapse_lines) ---
  contrasts_section <- if (!is.null(contrasts_text)) {
    c("EXPERIMENTAL CONTRASTS:", contrasts_text)
  }
  module_cors_section <- if (!is.null(module_cors_text) && length(module_cors_text) > 0) {
    c("MODULE-MODULE EIGENGENE CORRELATIONS (|r| >= 0.70):", module_cors_text)
  }

  text <- omicsai::collapse_lines(
    paste(lines, collapse = "\n"),
    contrasts_section,
    module_cors_section,
    sep = "\n\n"
  )

  list(
    text = text,
    data = list(
      experiment = experiment,
      organism = organism,
      n_samples = n_samples,
      n_features = n_features,
      n_wgcna_features = n_wgcna_features,
      modules = module_data
    )
  )
}


#' Build prompt parameters for a single WGCNA module summary
#'
#' Calls extract_module_data() for one module and formats the result into
#' the template parameters that the summary prompt expects.
#'
#' @param wgcna WGCNA results object
#' @param module Character; module name (e.g., "MEblue")
#' @param pgx PGX object
#'
#' @return Named list with template parameters:
#'   module, phenotypes, experiment, genesets, keygenes_section, module_stats
wgcna_build_summary_params <- function(wgcna, module, pgx) {
  annot <- wgcna$annot
  md <- extract_module_data(wgcna, module, pgx)
  M <- tryCatch(playbase:::wgcna.get_modTraits(wgcna), error = function(e) NULL)

  # --- Genesets ---
  ss <- ""
  gse <- md$gse
  if (is.data.frame(gse) && nrow(gse) > 0 &&
    all(c("geneset", "score", "q.value", "overlap") %in% colnames(gse))) {
    gse <- gse[order(gse$q.value, gse$score, na.last = TRUE), , drop = FALSE]
    sig <- gse[!is.na(gse$q.value) & gse$q.value < 0.05, , drop = FALSE]
    top_gse <- head(sig, 8)
    if (nrow(top_gse) > 0) {
      gset_df <- data.frame(
        Pathway = top_gse$geneset,
        Score = top_gse$score,
        `Q-value` = top_gse$q.value,
        Overlap = top_gse$overlap,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      gset_table <- omicsai::omicsai_format_mdtable(gset_df, formatters = list(
        Score = function(x) omicsai::omicsai_format_num(x, 2),
        `Q-value` = omicsai::omicsai_format_pvalue
      ))

      score_vals <- gse$score[!is.na(gse$score)]
      score_range <- if (length(score_vals) > 0) {
        paste0(
          omicsai::omicsai_format_num(min(score_vals), 2), "-",
          omicsai::omicsai_format_num(max(score_vals), 2),
          " (median: ", omicsai::omicsai_format_num(median(score_vals), 2), ")"
        )
      } else {
        "NA"
      }
      n_sig <- sum(gse$q.value < 0.05, na.rm = TRUE)

      ss <- omicsai::collapse_lines(
        "**Top Enriched Pathways (q < 0.05):**",
        gset_table,
        paste0(
          "**Score range:** ", score_range, "  \n",
          "**Total significant pathways:** ", n_sig, " of ", nrow(gse), " tested"
        ),
        sep = "\n\n"
      )
    }
  }
  if (ss == "" && length(md$fallback_sets) > 0) {
    ss <- paste0("- ", md$fallback_sets, collapse = "\n")
  }

  # --- Key genes section (multi-trait approach) ---
  # Only phenotype-derived traits for summary mode. When phenotypes is "None",
  # traits is empty and we fall through to the gene-list fallback.
  # extract_module_data also computes stats for top_trait as a fallback for
  # report mode — filter those out here.
  keygenes_section <- ""
  trait_stats <- md$trait_stats[intersect(names(md$trait_stats), md$traits)]

  if (length(trait_stats) > 0) {
    base_stats <- trait_stats[[1]]
    if (!"moduleMembership" %in% colnames(base_stats)) {
      base_stats$moduleMembership <- NA_real_
    }
    if (!"centrality" %in% colnames(base_stats)) {
      base_stats$centrality <- NA_real_
    }

    if ("moduleMembership" %in% colnames(base_stats)) {
      base_stats <- base_stats[order(-abs(base_stats$moduleMembership)), , drop = FALSE]
    } else if ("score" %in% colnames(base_stats)) {
      base_stats <- base_stats[order(-base_stats$score), , drop = FALSE]
    }

    top_genes <- head(base_stats, 8)
    features <- top_genes$feature
    symbols <- resolve_symbols(features, annot)

    # Table 1: trait-independent network metrics
    network_df <- data.frame(
      Gene = symbols,
      MM = top_genes$moduleMembership,
      Centrality = top_genes$centrality,
      stringsAsFactors = FALSE
    )
    network_table <- omicsai::collapse_lines(
      "### Trait-independent metrics (network structure)",
      omicsai::omicsai_format_mdtable(network_df, formatters = list(
        MM = function(x) omicsai::omicsai_format_num(x, 2),
        Centrality = function(x) omicsai::omicsai_format_num(x, 2)
      )),
      sep = "\n\n"
    )

    # Tables 2..N: per-trait gene-trait associations
    # Pre-compute module row index for M (constant across traits)
    mod_row_idx <- if (!is.null(M)) {
      which(rownames(M) %in% c(module, paste0("ME", module)))
    } else {
      integer(0)
    }

    trait_tables <- vapply(names(trait_stats), function(tr) {
      gs <- trait_stats[[tr]]

      mod_cor <- NA_real_
      if (length(mod_row_idx) > 0) {
        ti <- which(colnames(M) == tr)
        if (length(ti) > 0) mod_cor <- M[mod_row_idx[1], ti[1]]
      }
      cor_label <- if (!is.na(mod_cor)) {
        paste0(" (module r = ", omicsai::omicsai_format_num(mod_cor, 2), ")")
      } else {
        ""
      }

      ts_vals <- vapply(features, function(f) {
        idx <- match(f, gs$feature)
        if (!is.na(idx) && "traitSignificance" %in% colnames(gs)) gs$traitSignificance[idx] else NA_real_
      }, numeric(1))
      fc_vals <- vapply(features, function(f) {
        idx <- match(f, gs$feature)
        if (!is.na(idx) && "foldChange" %in% colnames(gs)) gs$foldChange[idx] else NA_real_
      }, numeric(1))

      trait_df <- data.frame(
        Gene = symbols, TS = ts_vals, logFC = fc_vals,
        stringsAsFactors = FALSE
      )
      omicsai::collapse_lines(
        paste0("### Gene-trait associations: ", tr, cor_label),
        omicsai::omicsai_format_mdtable(trait_df, formatters = list(
          TS = function(x) omicsai::omicsai_format_num(x, 2),
          logFC = function(x) omicsai::omicsai_format_num(x, 1)
        )),
        sep = "\n\n"
      )
    }, character(1))

    keygenes_section <- omicsai::collapse_lines(
      paste0(
        "**Hub Genes (ranked by module membership, top ",
        nrow(top_genes), " of ",
        if (!is.na(md$size)) md$size else nrow(base_stats),
        " module genes):**"
      ),
      network_table,
      paste(trait_tables, collapse = "\n\n"),
      paste(
        "**MM:** Module Membership (correlation with eigengene) — trait-independent  ",
        "**Centrality:** Intramodular connectivity — trait-independent  ",
        "**TS:** Trait Significance (gene-trait correlation) — specific to each trait  ",
        "**logFC:** Log2 fold change — specific to each trait",
        sep = "\n"
      ),
      sep = "\n\n"
    )
  } else if (length(md$fallback_genes) > 0) {
    genes <- head(md$fallback_genes, 15)
    keygenes_section <- paste0(
      "The following hub genes show high intramodular connectivity:\n\n",
      paste(genes, collapse = ", ")
    )
  }

  # --- Module statistics ---
  stat_lines <- c("**Module Statistics:**")
  if (!is.na(md$size)) {
    stat_lines <- c(stat_lines, paste0("- **Size:** ", md$size, " genes"))
  }

  if (length(trait_stats) > 0) {
    mod_row_idx_stats <- if (!is.null(M)) {
      which(rownames(M) %in% c(module, paste0("ME", module)))
    } else {
      integer(0)
    }
    for (tr in names(trait_stats)) {
      if (length(mod_row_idx_stats) > 0) {
        ti <- which(colnames(M) == tr)
        if (length(ti) > 0) {
          trait_cor <- M[mod_row_idx_stats[1], ti[1]]
          if (!is.na(trait_cor)) {
            stat_lines <- c(stat_lines, paste0(
              "- **Trait correlation with ", tr, ":** ", omicsai::omicsai_format_num(trait_cor, 2)
            ))
          }
        }
      }
      gs <- trait_stats[[tr]]
      if ("foldChange" %in% colnames(gs)) {
        mean_fc <- mean(gs$foldChange, na.rm = TRUE)
        if (!is.na(mean_fc)) {
          fc_dir <- ifelse(mean_fc > 0, "upregulated", "downregulated")
          stat_lines <- c(stat_lines, paste0(
            "- **Mean expression change (", tr, "):** ",
            omicsai::omicsai_format_num(abs(mean_fc), 1), "-fold ", fc_dir
          ))
        }
      }
    }
  }

  module_stats <- paste(stat_lines, collapse = "\n")

  list(
    module = module,
    phenotypes = md$phenotypes,
    experiment = wgcna$experiment %||% "",
    genesets = ss,
    keygenes_section = keygenes_section,
    module_stats = module_stats
  )
}


#' Rank WGCNA modules into signal tiers
#'
#' Deterministic classification of modules into strong / moderate / weak
#' tiers based on enrichment counts, trait correlations, and module size.
#'
#' @param wgcna WGCNA results object
#'
#' @return A list with elements: strong, moderate, weak, artifact,
#'   artifact_flags, order, grey_size, classification
wgcna_rank_modules <- function(wgcna) {
  M <- playbase:::wgcna.get_modTraits(wgcna)

  all_modules <- names(wgcna$me.genes)
  modules <- setdiff(all_modules, "MEgrey")

  grey_size <- length(wgcna$me.genes[["MEgrey"]])

  classification <- list()
  for (mod in modules) {
    # Count significant enrichments
    gse <- wgcna$gse[[mod]]
    if (!is.null(gse) && is.data.frame(gse) && "q.value" %in% colnames(gse)) {
      n_sig <- sum(gse$q.value < 0.05, na.rm = TRUE)
      n_total <- nrow(gse)
    } else {
      n_sig <- 0L
      n_total <- 0L
    }

    # Max absolute trait correlation
    mod_row <- mod
    if (mod_row %in% rownames(M)) {
      max_r <- max(abs(M[mod_row, ]), na.rm = TRUE)
    } else {
      max_r <- 0
    }

    # Module size
    size <- length(wgcna$me.genes[[mod]])

    # Classify tier
    if (n_sig >= 10 && max_r > 0.7 && size >= 50) {
      tier <- "strong"
    } else if (n_sig >= 1) {
      tier <- "moderate"
    } else {
      tier <- "weak"
    }

    classification[[mod]] <- list(
      tier = tier,
      n_sig = n_sig,
      max_r = max_r,
      size = size,
      n_total = n_total
    )
  }

  # Artifact detection
  artifact_flags <- list()
  for (mod in modules) {
    af <- classify_artifact(wgcna, mod)
    if (af$confidence != "none") {
      artifact_flags[[mod]] <- af
    }
  }

  artifact <- names(Filter(function(x) x$confidence == "probable", artifact_flags))

  # Build tier vectors
  strong <- names(Filter(function(x) x$tier == "strong", classification))
  moderate <- names(Filter(function(x) x$tier == "moderate", classification))
  weak <- names(Filter(function(x) x$tier == "weak", classification))

  # Order: strong first, then moderate, then weak
  order <- c(strong, moderate, weak)

  list(
    strong = strong,
    moderate = moderate,
    weak = weak,
    artifact = artifact,
    artifact_flags = artifact_flags,
    order = order,
    grey_size = as.integer(grey_size),
    classification = classification
  )
}


#' Build deterministic methods section for WGCNA report
#'
#' Reads the methods template and fills in WGCNA/PGX parameters using
#' omicsai::omicsai_substitute_template().
#'
#' @param wgcna WGCNA results object
#' @param pgx PGX object
#'
#' @return Character string with the rendered methods section
wgcna_build_methods <- function(wgcna, pgx) {
  template_path <- file.path(BOARD_PROMPTS_DIR, "wgcna_report_methods.md")
  template <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

  # Extract values (same logic as prompt_optim/pipeline/build_methods.R)
  n_wgcna <- length(unlist(wgcna$me.genes))
  if (is.null(n_wgcna) || n_wgcna == 0) {
    n_wgcna <- tryCatch(ncol(wgcna$datExpr), error = function(e) NA)
  }

  feature_type <- "features"
  if (!is.null(pgx$datatype)) {
    if (grepl("prot", pgx$datatype, ignore.case = TRUE)) {
      feature_type <- "proteins"
    } else if (grepl("rna|transcript|gene", pgx$datatype, ignore.case = TRUE)) {
      feature_type <- "genes"
    }
  }

  n_samples <- tryCatch(nrow(wgcna$datExpr), error = function(e) {
    tryCatch(nrow(pgx$samples), error = function(e2) NA)
  })

  network_type <- wgcna$net$networkType %||% wgcna$networkType %||% "signed"
  power <- wgcna$power %||% wgcna$net$power %||% "NA"
  min_mod_size <- wgcna$minModSize %||% "20"
  merge_cut_height <- wgcna$mergeCutHeight %||% "0.15"
  min_kme <- wgcna$minKME %||% "0.3"

  modules_no_grey <- setdiff(names(wgcna$me.genes), c("grey", "MEgrey"))
  n_modules <- length(modules_no_grey)
  grey_size <- length(wgcna$me.genes[["MEgrey"]]) + length(wgcna$me.genes[["grey"]])

  n_genesets <- tryCatch({
    first_mod <- modules_no_grey[1]
    nrow(wgcna$gse[[first_mod]])
  }, error = function(e) "NA")

  params <- list(
    n_features_wgcna = as.character(n_wgcna),
    feature_type = feature_type,
    n_samples = as.character(n_samples),
    network_type = as.character(network_type),
    power = as.character(power),
    min_mod_size = as.character(min_mod_size),
    merge_cut_height = as.character(merge_cut_height),
    min_kme = as.character(min_kme),
    n_modules = as.character(n_modules),
    grey_size = as.character(grey_size),
    n_genesets_tested = as.character(n_genesets),
    date = format(Sys.Date(), "%Y-%m-%d")
  )

  omicsai::omicsai_substitute_template(template, params)
}
