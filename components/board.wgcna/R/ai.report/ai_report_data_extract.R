## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.

# =============================================================================
# WGCNA AI Report — Extraction Helpers
# =============================================================================
# Private helpers used by the public functions in ai_report_data.R.
# No Shiny dependencies. No LLM calls.
#
# Functions:
#   extract_module_data(wgcna, module, pgx, ...)  — shared per-module extraction core
#   resolve_symbols(features, annot)               — feature IDs to gene symbols
#   resolve_functions(features, annot)             — feature IDs to descriptions
#   select_top_enrichment(sig_gse, n)              — tiered enrichment term selection
#   classify_artifact(wgcna, module)               — serum/hemoglobin contamination check

# -----------------------------------------------------------------------------
# Symbol and annotation resolution
# -----------------------------------------------------------------------------

#' Resolve feature IDs to gene symbols
#'
#' @param features Character vector of feature IDs
#' @param annot Annotation data frame (NULL-safe)
#' @return Character vector of symbols (same length as features)
resolve_symbols <- function(features, annot) {
  if (is.null(annot)) return(features)
  symbols <- tryCatch(
    playbase::probe2symbol(features, annot, "symbol"),
    error = function(e) features
  )
  ifelse(is.na(symbols) | symbols == "", features, symbols)
}

#' Look up function descriptions for features
#'
#' @param features Character vector of feature IDs
#' @param annot Annotation data frame (NULL-safe)
#' @param max_chars Integer; truncation limit
#' @return Character vector of descriptions
resolve_functions <- function(features, annot, max_chars = 60L) {
  funcs <- rep("", length(features))
  if (is.null(annot)) return(funcs)
  func_col <- intersect(c("gene_title", "gene_name", "description"), colnames(annot))
  if (length(func_col) == 0) return(funcs)
  idx <- match(features, rownames(annot))
  valid <- !is.na(idx)
  funcs[valid] <- as.character(annot[idx[valid], func_col[1]])
  substr(funcs, 1, max_chars)
}

# -----------------------------------------------------------------------------
# Enrichment selection
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

# -----------------------------------------------------------------------------
# Artifact detection
# -----------------------------------------------------------------------------

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

  symbols <- resolve_symbols(hub_features, wgcna$annot)

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
# Shared extraction core
# -----------------------------------------------------------------------------

#' Extract all data for a single WGCNA module
#'
#' Shared extraction core used by both wgcna_build_report_tables (report mode)
#' and wgcna_build_summary_params (summary mode). Extracts phenotypes,
#' enrichment, hub gene stats, trait correlations, and eigengene profiles.
#'
#' @param wgcna WGCNA results object
#' @param module Character; module name (e.g. "MEblue")
#' @param pgx PGX object
#' @param top Pre-computed result from wgcna.getTopGenesAndSets (optional;
#'   computed if NULL). Pass this when calling in a loop to avoid recomputation.
#' @param M Pre-computed module-trait correlation matrix (optional;
#'   computed if NULL). Pass this when calling in a loop to avoid recomputation.
#' @param max_traits Integer; max correlated traits for gene stats (default 3)
#'
#' @return A list with raw extracted data for the module.
#'
#' @note Multi-omics WGCNA (consensusWGCNA, multiwgcna) will have their own
#'   board-specific extraction; this function handles standard WGCNA only.
extract_module_data <- function(wgcna, module, pgx,
                                top = NULL, M = NULL,
                                max_traits = 3L) {
  annot <- wgcna$annot

  # --- Phenotypes via wgcna.getTopGenesAndSets ---
  if (is.null(top)) {
    top <- tryCatch(
      playbase::wgcna.getTopGenesAndSets(wgcna, annot = annot, ntop = 40,
                                         level = "gene", rename = "gene_title"),
      error = function(e) list(pheno = list(), genes = list(), sets = list())
    )
  }

  phenotypes <- "None"
  if (module %in% names(top$pheno)) {
    phenotypes <- paste(top$pheno[[module]], collapse = ", ")
  }
  fallback_genes <- top$genes[[module]] %||% character(0)
  fallback_sets <- top$sets[[module]] %||% character(0)

  # --- Module size ---
  size <- length(wgcna$me.genes[[module]])

  # --- Module-trait correlation matrix ---
  if (is.null(M)) {
    M <- tryCatch(playbase:::wgcna.get_modTraits(wgcna), error = function(e) NULL)
  }

  top_trait <- ""
  top_r <- NA_real_
  if (!is.null(M) && module %in% rownames(M)) {
    cors <- M[module, ]
    top_idx <- which.max(abs(cors))
    if (length(top_idx) > 0) {
      top_trait <- names(cors)[top_idx]
      top_r <- cors[top_idx]
    }
  }

  # --- Correlated traits (from phenotypes) ---
  traits <- character(0)
  if (phenotypes != "None" && nzchar(phenotypes)) {
    traits <- trimws(strsplit(phenotypes, ",\\s*")[[1]])
    traits <- traits[nzchar(traits) & traits != "None"]
    traits <- head(traits, max_traits)
  }

  # --- Eigengene profile ---
  eigengene_profile <- NULL
  peak_condition <- ""
  datME <- wgcna$datME %||% wgcna$net$MEs
  groups <- if (!is.null(pgx$samples$group)) pgx$samples$group else NULL
  if (!is.null(datME) && module %in% colnames(datME) && !is.null(groups)) {
    eigengene <- datME[, module]
    group_means <- tapply(eigengene, groups, mean)
    eigengene_profile <- group_means
    peak_condition <- names(which.max(abs(group_means)))
  }

  # --- Enrichment data ---
  n_sig <- 0L
  n_total <- 0L
  gse <- wgcna$gse[[module]]
  if (is.data.frame(gse) && nrow(gse) > 0 && "q.value" %in% colnames(gse)) {
    n_sig <- sum(gse$q.value < 0.05, na.rm = TRUE)
    n_total <- nrow(gse)
  }

  # --- Gene stats per trait ---
  trait_stats <- list()
  # Use phenotype-derived traits if available, otherwise fall back to top trait
  stat_traits <- if (length(traits) > 0) traits else if (nzchar(top_trait)) top_trait
  for (tr in stat_traits) {
    gs <- tryCatch(
      playbase::wgcna.getGeneStats(wgcna, module = module, trait = tr, plot = FALSE),
      error = function(e) NULL
    )
    if (is.data.frame(gs) && nrow(gs) > 0) {
      trait_stats[[tr]] <- gs
    }
  }

  list(
    size = size,
    phenotypes = phenotypes,
    traits = traits,
    fallback_genes = fallback_genes,
    fallback_sets = fallback_sets,
    eigengene_profile = eigengene_profile,
    peak_condition = peak_condition,
    top_trait = top_trait,
    top_r = top_r,
    n_sig = n_sig,
    n_total = n_total,
    gse = gse,
    trait_stats = trait_stats
  )
}
