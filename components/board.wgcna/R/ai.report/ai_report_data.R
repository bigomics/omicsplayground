## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.

# =============================================================================
# WGCNA AI Report — Data Extraction and Classification
# =============================================================================
# Deterministic functions that build structured data for the AI report prompt.
# No Shiny dependencies. No LLM calls.
#
# Public functions:
#   wgcna_build_report_tables(wgcna, pgx, ...)  — structured data tables
#   wgcna_rank_modules(wgcna)                    — signal tier classification
#   wgcna_build_methods(wgcna, pgx)              — deterministic methods section

BOARD_PROMPTS_DIR <- file.path(OPG, "components/board.wgcna/prompts")

# -----------------------------------------------------------------------------
# Local helpers
# -----------------------------------------------------------------------------

#' Select top enrichment terms with tiered preference
#'
#' @param sig_gse data.frame of significant enrichment terms, pre-sorted by q.value
#' @param n Integer; max terms to return
#' @return data.frame with up to n rows
select_top_enrichment <- function(sig_gse, n) {
  # Deduplicate GO terms: extract GO accession, keep first occurrence
  go_id <- ifelse(
    grepl("GO_\\d{5,}", sig_gse$geneset),
    sub(".*?(GO_\\d{5,}).*", "\\1", sig_gse$geneset),
    NA_character_
  )
  dup <- duplicated(go_id) & !is.na(go_id)
  sig_gse <- sig_gse[!dup, , drop = FALSE]

  prefix <- sub(":.*", "", sig_gse$geneset)

  # Tier 1: biological process / pathway databases
  tier1 <- c("HALLMARK", "C5", "GO_BP", "GOBP", "GO_CC", "GOCC", "GO_MF", "GOMF",
             "PATHWAY", "PATHWAY_REACTOME", "PATHWAY_KEGG", "PATHWAY_BIOPLANET",
             "PATHWAY_WIKI")
  # Tier 2: regulatory mechanisms
  tier2 <- c("TF_ARCHS4", "KINASE_ARCHS4")

  idx1 <- which(prefix %in% tier1)
  idx2 <- which(prefix %in% tier2)
  idx_rest <- setdiff(seq_len(nrow(sig_gse)), c(idx1, idx2))

  # Fill from tier 1 first, then tier 2, then remainder
  selected <- head(c(idx1, idx2, idx_rest), n)
  selected <- sort(selected)  # preserve q-value ordering within selection

  sig_gse[selected, , drop = FALSE]
}

#' Classify artifact contamination in a WGCNA module
#'
#' @param wgcna WGCNA results object
#' @param module Character; module name (e.g. "MEblue")
#' @return A list with confidence ("probable", "possible", "none") and reason.
classify_artifact <- function(wgcna, module) {
  SERUM_MARKERS <- c(
    "ALB", "C3", "C4A", "C4B", "C5", "FGA", "FGB", "FGG",
    "SERPINA1", "APOA1", "APOB", "HP", "HPX", "TF", "AHSG",
    "FN1", "A2M", "CFH", "CFB", "ORM1", "GC"
  )
  HEMOGLOBIN_MARKERS <- c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2")
  ALL_MARKERS <- c(SERUM_MARKERS, HEMOGLOBIN_MARKERS)

  # Get hub genes via trait-based gene stats
  gs <- NULL
  M <- tryCatch(playbase:::wgcna.get_modTraits(wgcna), error = function(e) NULL)
  if (!is.null(M)) {
    mod_row <- module
    if (mod_row %in% rownames(M)) {
      cors <- M[mod_row, ]
      top_trait <- names(which.max(abs(cors)))
      gs <- tryCatch(
        playbase::wgcna.getGeneStats(
          wgcna,
          module = module,
          trait = top_trait,
          plot = FALSE
        ),
        error = function(e) NULL
      )
    }
  }

  # Extract top 10 hub genes by module membership
  hub_features <- NULL
  if (is.data.frame(gs) && nrow(gs) > 0) {
    if ("moduleMembership" %in% colnames(gs)) {
      gs <- gs[order(-abs(gs$moduleMembership)), , drop = FALSE]
    }
    hub_features <- head(gs$feature, 10)
  }

  # Fallback: use module gene list directly
  if (is.null(hub_features) || length(hub_features) == 0) {
    hub_features <- head(wgcna$me.genes[[module]], 10)
  }

  if (is.null(hub_features) || length(hub_features) == 0) {
    return(list(confidence = "none", reason = NULL))
  }

  # Resolve symbols
  symbols <- hub_features
  if (!is.null(wgcna$annot)) {
    symbols <- tryCatch(
      playbase::probe2symbol(hub_features, wgcna$annot, "symbol"),
      error = function(e) hub_features
    )
    symbols <- ifelse(is.na(symbols) | symbols == "", hub_features, symbols)
  }

  matches <- intersect(toupper(symbols), ALL_MARKERS)
  n_matches <- length(matches)

  if (n_matches >= 3) {
    list(
      confidence = "probable",
      reason = paste("serum proteins in hub genes:", paste(matches, collapse = ", "))
    )
  } else if (n_matches >= 1) {
    list(
      confidence = "possible",
      reason = paste("serum proteins in hub genes:", paste(matches, collapse = ", "))
    )
  } else {
    list(confidence = "none", reason = NULL)
  }
}

# -----------------------------------------------------------------------------
# Public functions
# -----------------------------------------------------------------------------

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

  # --- Module-trait correlation matrix ---
  M <- tryCatch(playbase:::wgcna.get_modTraits(wgcna), error = function(e) NULL)

  # --- Eigengene data ---
  datME <- wgcna$datME %||% wgcna$net$MEs
  groups <- if (!is.null(pgx$samples$group)) pgx$samples$group else NULL

  # --- Build per-module data ---
  module_data <- list()
  overview_rows <- list()

  for (mod in non_grey) {
    size <- length(wgcna$me.genes[[mod]])

    # Eigengene profile
    eigengene_profile <- NULL
    peak_condition <- ""
    me_col <- mod
    if (!is.null(datME) && me_col %in% colnames(datME) && !is.null(groups)) {
      eigengene <- datME[, me_col]
      group_means <- tapply(eigengene, groups, mean)
      eigengene_profile <- group_means
      peak_condition <- names(which.max(abs(group_means)))
    }

    # Top trait from M
    top_trait <- ""
    top_r <- NA_real_
    mod_row <- mod
    if (!is.null(M) && mod_row %in% rownames(M)) {
      cors <- M[mod_row, ]
      top_idx <- which.max(abs(cors))
      if (length(top_idx) > 0) {
        top_trait <- names(cors)[top_idx]
        top_r <- cors[top_idx]
      }
    }

    # Enrichment data — tiered selection preferring interpretable terms
    n_sig <- 0L
    n_total <- 0L
    top_enrich <- NULL
    gse <- wgcna$gse[[mod]]
    if (is.data.frame(gse) && nrow(gse) > 0 && "q.value" %in% colnames(gse)) {
      n_sig <- sum(gse$q.value < 0.05, na.rm = TRUE)
      n_total <- nrow(gse)
      sig_gse <- gse[!is.na(gse$q.value) & gse$q.value < 0.05, , drop = FALSE]
      if (nrow(sig_gse) > 0) {
        sig_gse <- sig_gse[order(sig_gse$q.value), , drop = FALSE]
        top_enrich <- select_top_enrichment(sig_gse, ntop_enrichment)
      }
    }

    # Hub genes via wgcna_build_ai_params for trait info, then gene stats
    ai_params <- tryCatch(
      wgcna_build_ai_params(wgcna, mod, annot = wgcna$annot),
      error = function(e) NULL
    )

    hub_genes_df <- NULL
    trait_for_genes <- top_trait
    if (!is.null(ai_params) && ai_params$phenotypes != "None") {
      traits <- trimws(strsplit(ai_params$phenotypes, ",\\s*")[[1]])
      traits <- traits[nzchar(traits) & traits != "None"]
      if (length(traits) > 0) trait_for_genes <- traits[1]
    }

    if (nzchar(trait_for_genes)) {
      gs <- tryCatch(
        playbase::wgcna.getGeneStats(wgcna, module = mod,
                                     trait = trait_for_genes, plot = FALSE),
        error = function(e) NULL
      )
      if (is.data.frame(gs) && nrow(gs) > 0) {
        if ("moduleMembership" %in% colnames(gs)) {
          gs <- gs[order(-abs(gs$moduleMembership)), , drop = FALSE]
        }
        top_gs <- head(gs, ntop_genes)

        # Resolve symbols
        symbols <- top_gs$feature
        if (!is.null(wgcna$annot)) {
          symbols <- tryCatch(
            playbase::probe2symbol(top_gs$feature, wgcna$annot, "symbol"),
            error = function(e) top_gs$feature
          )
          symbols <- ifelse(is.na(symbols) | symbols == "", top_gs$feature, symbols)
        }

        # Get function descriptions
        funcs <- rep("", nrow(top_gs))
        if (!is.null(wgcna$annot)) {
          func_col <- intersect(c("gene_title", "gene_name", "description"),
                                colnames(wgcna$annot))
          if (length(func_col) > 0) {
            idx <- match(top_gs$feature, rownames(wgcna$annot))
            valid <- !is.na(idx)
            funcs[valid] <- as.character(wgcna$annot[idx[valid], func_col[1]])
            funcs <- substr(funcs, 1, 60)
          }
        }

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
    }

    module_data[[mod]] <- list(
      size = size,
      eigengene_profile = eigengene_profile,
      peak_condition = peak_condition,
      top_trait = top_trait,
      top_r = top_r,
      n_sig = n_sig,
      n_total = n_total,
      top_enrichment = top_enrich,
      hub_genes = hub_genes_df,
      trait_for_genes = trait_for_genes
    )

    overview_rows[[mod]] <- list(
      module = mod,
      size = size,
      peak_condition = peak_condition,
      top_trait = top_trait,
      top_r = top_r,
      n_sig = n_sig,
      n_total = n_total
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
      md <- module_data[[mod]]
      if (!is.null(md$top_enrichment) && nrow(md$top_enrichment) > 0) {
        te <- md$top_enrichment
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
    md <- module_data[[mod]]
    lines <- c(lines, sprintf("### %s (%d genes)", mod, md$size))

    # Eigengene profile
    if (!is.null(md$eigengene_profile)) {
      profile_str <- paste(
        sprintf("%s=%s", names(md$eigengene_profile),
                omicsai::omicsai_format_num(md$eigengene_profile, 2)),
        collapse = ", "
      )
      lines <- c(lines, paste0("Eigengene profile: ", profile_str))
    }

    # Top trait
    if (nzchar(md$top_trait) && !is.na(md$top_r)) {
      lines <- c(lines, sprintf("Top trait: %s (r=%+.2f)", md$top_trait, md$top_r))
    }
    lines <- c(lines, "")

    # Enrichment
    if (!is.null(md$top_enrichment) && nrow(md$top_enrichment) > 0) {
      te <- md$top_enrichment
      enrich_df <- data.frame(
        `#` = seq_len(nrow(te)),
        Term = te$geneset,
        `q-value` = te$q.value,
        Source = if ("source" %in% colnames(te)) as.character(te$source) else "",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      lines <- c(lines,
        sprintf("Top enrichment (%d significant of %d):", md$n_sig, md$n_total),
        omicsai::omicsai_format_mdtable(enrich_df, formatters = list(
          `q-value` = omicsai::omicsai_format_pvalue
        ))
      )
    } else {
      lines <- c(lines, "No significant enrichment (all q > 0.05)")
    }
    lines <- c(lines, "")

    # Hub genes
    if (!is.null(md$hub_genes) && nrow(md$hub_genes) > 0) {
      hg <- md$hub_genes
      logfc_label <- if (nzchar(md$trait_for_genes)) {
        paste0("logFC (vs ", md$trait_for_genes, ")")
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
    n_genesets_tested = as.character(n_genesets)
  )

  omicsai::omicsai_substitute_template(template, params)
}
