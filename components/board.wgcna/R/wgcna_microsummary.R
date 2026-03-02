##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# WGCNA Microsummary — board-specific static texts and param extractor
# =============================================================================
# Consumed by MicrosummaryServer() in wgcna_server.R.
# Each board that adopts microsummaries creates its own <board>_microsummary.R
# with the same two exports: <BOARD>_STATIC_TEXTS + <board>_microsummary_params().

#' Static descriptions for each WGCNA tab (originally inline in wgcna_ui.R)
WGCNA_STATIC_TEXTS <- list(
  "WGCNA" = '<b>Module detection.</b> <b>(a)</b> Modules are detected using the dynamic branch cutting approach. <b>(b)</b> Scale independence and mean connectivity plots to determine the soft threshold. <b>(c)</b> Topological overlap matrix visualized as heatmap. <b>(d)</b> Dimensionality reduction map of features colored by module. <b>(e)</b> Size of WGCNA modules.',
  "Eigengenes" = '<b>Eigengene analysis.</b> The module eigengene of a given module is defined as the first principal component of the standardized expression profiles. <b>(a)</b> Module-trait correlation identifies modules that are significantly associated with the measured traits. <b>(b)</b> Clustering of eigengenes. <b>(c)</b> Clustering of trait vectors. <b>(d)</b> Correlation of eigengene and traits as heatmap, <b>(e)</b> as dendrogram and <b>(f)</b> as graph.',
  "Modules" = '<b>Module analysis.</b>  <b>(a)</b> Correlation of module eigengene with traits. <b>(b)</b> Circle network of top hub genes. <b>(c)</b> Table of importance score to identify \'driver genes\' of the module. <b>(d)</b> Plot of gene significance parameters.',
  "Enrichment" = '<b>Module Enrichment.</b> <b>(a)</b> Enrichment heatmap of top most enriched genesets in module. <b>(b)</b> Expression heatmap of genes in selected geneset. <b>(c)</b> Functional enrichment of the module calculated using Fisher\'s exact test. <b>(d)</b> Top enriched genesets in module.'
)

#' Extract light parameters for WGCNA microsummary
#'
#' Returns a named list with dataset_name and tab_context for each WGCNA tab.
#' Uses only lightweight extraction (no wgcna.getGeneStats).
#'
#' @param tab Character; tab name ("WGCNA", "Eigengenes", "Modules", "Enrichment")
#' @param wgcna WGCNA results object
#' @param pgx PGX object
#' @param selected_module Character; currently selected module (for Modules/Enrichment tabs)
#' @return Named list with dataset_name and tab_context, or NULL if data missing
wgcna_microsummary_params <- function(tab, wgcna, pgx, selected_module = NULL) {
  if (is.null(wgcna) || is.null(wgcna$me.genes)) return(NULL)

  dataset_name <- wgcna$experiment %||% pgx$name %||% pgx$description %||% "dataset"

  context <- switch(tab,
    "WGCNA" = .wgcna_micro_tab_wgcna(wgcna, pgx),
    "Eigengenes" = .wgcna_micro_tab_eigengenes(wgcna, pgx),
    "Modules" = .wgcna_micro_tab_modules(wgcna, pgx, selected_module),
    "Enrichment" = .wgcna_micro_tab_enrichment(wgcna, pgx, selected_module),
    NULL
  )

  if (is.null(context)) return(NULL)

  list(
    dataset_name = dataset_name,
    tab_context = context
  )
}

# -- Per-tab context builders --------------------------------------------------

.wgcna_micro_tab_wgcna <- function(wgcna, pgx) {
  all_modules <- names(wgcna$me.genes)
  non_grey <- setdiff(all_modules, "MEgrey")
  n_modules <- length(non_grey)
  total_genes <- length(unlist(wgcna$me.genes))
  power <- wgcna$power %||% wgcna$net$power %||% "auto"
  n_grey <- length(wgcna$me.genes[["MEgrey"]])

  # Top 3 largest modules
  sizes <- lengths(wgcna$me.genes[non_grey])
  top3 <- head(sort(sizes, decreasing = TRUE), 3)
  top3_str <- paste(sprintf("%s (%d genes)", names(top3), top3), collapse = ", ")

  paste(
    sprintf("Module detection found %d modules from %d genes (soft threshold power=%s).", n_modules, total_genes, power),
    sprintf("Largest modules: %s.", top3_str),
    sprintf("%d genes were unassigned (grey module).", n_grey),
    sep = "\n"
  )
}

.wgcna_micro_tab_eigengenes <- function(wgcna, pgx) {
  M <- tryCatch(playbase:::wgcna.get_modTraits(wgcna), error = function(e) NULL)
  if (is.null(M)) return("Module-trait correlations are not available.")

  n_samples <- tryCatch(nrow(pgx$samples), error = function(e) NA_integer_)
  groups <- if (!is.null(pgx$samples$group)) {
    paste(unique(pgx$samples$group), collapse = ", ")
  } else {
    "unknown"
  }

  # Top 3 trait-correlated modules (by max abs correlation)
  max_cors <- apply(abs(M), 1, max, na.rm = TRUE)
  max_cors <- sort(max_cors, decreasing = TRUE)
  top3 <- head(max_cors, 3)
  top3_details <- vapply(names(top3), function(mod) {
    cors <- M[mod, ]
    best_idx <- which.max(abs(cors))
    sprintf("%s (r=%.2f with %s)", mod, cors[best_idx], names(cors)[best_idx])
  }, character(1))

  paste(
    sprintf("Eigengene analysis across %s samples in groups: %s.", n_samples, groups),
    sprintf("Top trait-correlated modules: %s.", paste(top3_details, collapse = "; ")),
    sep = "\n"
  )
}

.wgcna_micro_tab_modules <- function(wgcna, pgx, selected_module) {
  if (is.null(selected_module) || !selected_module %in% names(wgcna$me.genes)) {
    return(NULL)
  }

  mod_genes <- wgcna$me.genes[[selected_module]]
  mod_size <- length(mod_genes)

  # Top 5 hub genes (by position in module gene list, which is ordered by membership)
  annot <- wgcna$annot
  hub_genes <- head(mod_genes, 5)
  hub_symbols <- tryCatch({
    syms <- playbase::probe2symbol(hub_genes, annot, "symbol")
    ifelse(is.na(syms) | syms == "", hub_genes, syms)
  }, error = function(e) hub_genes)

  # Top enriched pathways from gse
  gse <- wgcna$gse[[selected_module]]
  top_pathways <- ""
  if (is.data.frame(gse) && nrow(gse) > 0 && "q.value" %in% colnames(gse)) {
    sig <- gse[!is.na(gse$q.value) & gse$q.value < 0.05, , drop = FALSE]
    sig <- sig[order(sig$q.value), , drop = FALSE]
    top3 <- head(sig, 3)
    if (nrow(top3) > 0) {
      top_pathways <- paste(sprintf("%s (q=%.1e)", top3$geneset, top3$q.value), collapse = "; ")
    }
  }

  lines <- c(
    sprintf("Module %s contains %d genes.", selected_module, mod_size),
    sprintf("Top hub genes: %s.", paste(hub_symbols, collapse = ", "))
  )
  if (nzchar(top_pathways)) {
    lines <- c(lines, sprintf("Top enriched pathways: %s.", top_pathways))
  }
  paste(lines, sep = "\n")
}

.wgcna_micro_tab_enrichment <- function(wgcna, pgx, selected_module) {
  if (is.null(selected_module) || !selected_module %in% names(wgcna$me.genes)) {
    return(NULL)
  }

  gse <- wgcna$gse[[selected_module]]
  if (!is.data.frame(gse) || nrow(gse) == 0 || !"q.value" %in% colnames(gse)) {
    return(sprintf("Module %s: no enrichment data available.", selected_module))
  }

  sig <- gse[!is.na(gse$q.value) & gse$q.value < 0.05, , drop = FALSE]
  sig <- sig[order(sig$q.value), , drop = FALSE]
  n_sig <- nrow(sig)
  n_total <- nrow(gse)

  top5 <- head(sig, 5)
  if (nrow(top5) > 0) {
    geneset_lines <- paste(
      sprintf("  - %s (q=%.1e)", top5$geneset, top5$q.value),
      collapse = "\n"
    )
  } else {
    geneset_lines <- "  - No significant enrichments found (q < 0.05)"
  }

  paste(
    sprintf("Module %s enrichment analysis: %d significant of %d tested genesets.", selected_module, n_sig, n_total),
    sprintf("Top enriched genesets:\n%s", geneset_lines),
    sep = "\n"
  )
}
