##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for clustering summary
#'
#' Extracts parameters from clustering results to populate the prompt template.
#'
#' @param pgx PGX object
#' @param top_matrix List from getTopMatrix() with mat, idx, grp, annot, etc.
#' @param clust_annot_cor Matrix from getClustAnnotCorrelation(); correlation of
#'   gene modules with annotation terms (rows = terms, cols = modules)
#' @param cluster_method Character; clustering layout method (e.g., "umap", "pca")
#' @param ntop Integer; number of top annotation terms per module (default 10)
#'
#' @return Named list with template parameters:
#'   experiment, cluster_method, n_modules, cluster_info, enriched_terms, top_genes
clustering_build_ai_params <- function(pgx,
                                       top_matrix,
                                       clust_annot_cor = NULL,
                                       cluster_method = "umap",
                                       ntop = 10) {
  # Extract experiment description
  experiment <- ""
  if (!is.null(pgx$name) && nzchar(pgx$name)) {
    experiment <- pgx$name
  } else if (!is.null(pgx$description) && nzchar(pgx$description)) {
    experiment <- pgx$description
  }

  # --------------------------------------------------------------------------
  # Gene module composition from top_matrix
  # --------------------------------------------------------------------------

  mat <- top_matrix$mat
  idx <- top_matrix$idx
  grp <- top_matrix$grp
  n_modules <- length(unique(idx))

  cluster_info <- ""
  if (!is.null(mat) && !is.null(idx)) {
    module_ids <- sort(unique(idx))
    module_rows <- list()
    for (mid in module_ids) {
      genes_in_mod <- rownames(mat)[idx == mid]
      n_genes <- length(genes_in_mod)
      # Show up to 10 representative gene names
      show_genes <- head(genes_in_mod, 10)
      gene_str <- paste(show_genes, collapse = ", ")
      if (n_genes > 10) gene_str <- paste0(gene_str, ", ...")
      module_rows[[mid]] <- sprintf("| %s | %s | %s |", mid, n_genes, gene_str)
    }

    cluster_info <- paste0(
      "**Gene modules (", n_modules, " modules, ", nrow(mat), " total genes, ",
      ncol(mat), " samples):**\n\n",
      "| Module | N genes | Representative genes |\n",
      "|--------|---------|---------------------|\n",
      paste(unlist(module_rows), collapse = "\n")
    )

    # Add sample group info if available
    if (!is.null(grp) && length(grp) > 0) {
      grp_table <- table(grp)
      grp_lines <- sprintf("- **%s:** %d samples", names(grp_table), as.integer(grp_table))
      cluster_info <- paste0(
        cluster_info, "\n\n",
        "**Sample groups:**\n",
        paste(grp_lines, collapse = "\n")
      )
    }
  }

  # --------------------------------------------------------------------------
  # Functional annotation of gene modules
  # --------------------------------------------------------------------------

  enriched_terms <- ""
  if (!is.null(clust_annot_cor) && nrow(clust_annot_cor) > 0 && ncol(clust_annot_cor) > 0) {
    rho <- clust_annot_cor
    module_sections <- list()

    for (j in seq_len(ncol(rho))) {
      mod_name <- colnames(rho)[j]
      # Get top annotation terms by correlation for this module
      vals <- rho[, j]
      top_idx <- head(order(-vals), ntop)
      top_terms <- rownames(rho)[top_idx]
      top_vals <- vals[top_idx]

      # Filter to positive correlations only
      keep <- top_vals > 0
      if (sum(keep) == 0) next
      top_terms <- top_terms[keep]
      top_vals <- top_vals[keep]

      term_table <- paste0(
        "**", mod_name, ":**\n\n",
        "| Annotation term | Correlation (R) |\n",
        "|----------------|----------------|\n",
        paste(
          sprintf("| %s | %s |", top_terms, omicsai::format_num(top_vals, 3)),
          collapse = "\n"
        )
      )
      module_sections <- c(module_sections, term_table)
    }

    if (length(module_sections) > 0) {
      enriched_terms <- paste(module_sections, collapse = "\n\n")
    }
  }

  if (enriched_terms == "") {
    enriched_terms <- "No functional annotation data available for the current clustering configuration."
  }

  # --------------------------------------------------------------------------
  # Top marker genes per module
  # --------------------------------------------------------------------------

  top_genes <- ""
  if (!is.null(mat) && !is.null(idx) && ncol(mat) > 1) {
    module_ids <- sort(unique(idx))
    gene_sections <- list()

    for (mid in module_ids) {
      mod_genes <- rownames(mat)[idx == mid]
      if (length(mod_genes) == 0) next

      # Compute mean expression and SD for each gene in this module
      mod_mat <- mat[mod_genes, , drop = FALSE]
      gene_sd <- apply(mod_mat, 1, sd, na.rm = TRUE)
      gene_mean <- rowMeans(mod_mat, na.rm = TRUE)

      # Rank by SD (most variable = most informative markers)
      ord <- order(-gene_sd)
      show_n <- min(15, length(ord))
      top_g <- mod_genes[ord[seq_len(show_n)]]

      gene_table <- paste0(
        "**", mid, "** (top ", show_n, " by variability):\n\n",
        "| Gene | Mean expr | SD |\n",
        "|------|-----------|----|\n",
        paste(
          sprintf(
            "| %s | %s | %s |",
            top_g,
            omicsai::format_num(gene_mean[top_g], 2),
            omicsai::format_num(gene_sd[top_g], 2)
          ),
          collapse = "\n"
        )
      )
      gene_sections <- c(gene_sections, gene_table)
    }

    if (length(gene_sections) > 0) {
      top_genes <- paste(gene_sections, collapse = "\n\n")
    }
  }

  if (top_genes == "") {
    top_genes <- "No marker gene data available for the current clustering configuration."
  }

  list(
    experiment = experiment,
    cluster_method = cluster_method,
    n_modules = as.character(n_modules),
    cluster_info = cluster_info,
    enriched_terms = enriched_terms,
    top_genes = top_genes
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: Clustering AI Summary
# -----------------------------------------------------------------------------

#' Clustering AI Summary Server
#'
#' Server logic for clustering AI-generated summaries.
#' Uses omicsai_summary_card_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param getTopMatrix Reactive returning top matrix list (mat, idx, grp, annot)
#' @param getClustAnnotCorrelation Reactive returning annotation correlation matrix
#' @param clustmethod Reactive returning current clustering method name
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_summary_card_server
clustering_ai_summary_server <- function(id,
                                         pgx,
                                         getTopMatrix,
                                         getClustAnnotCorrelation,
                                         clustmethod,
                                         session,
                                         cache = NULL,
                                         watermark = FALSE) {

  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    top_matrix <- getTopMatrix()
    shiny::req(top_matrix, top_matrix$mat)

    # Get annotation correlation (may be NULL if not yet computed)
    clust_annot_cor <- tryCatch(
      getClustAnnotCorrelation(),
      error = function(e) NULL
    )

    method <- clustmethod()
    shiny::req(method)

    clustering_build_ai_params(
      pgx = pgx,
      top_matrix = top_matrix,
      clust_annot_cor = clust_annot_cor,
      cluster_method = method
    )
  })

  # Load template from board.clustering prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.clustering/prompts/clustering_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omicsai::omicsai_load_template(board_template_path)
  })

  # Load clustering methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.clustering/prompts/Clustering_methods.md"
  )
  context_template <- omicsai::omicsai_load_template(context_template_path)

  # Build config reactive from user's selected LLM model
  # Includes board-specific context in system prompt
  config_reactive <- shiny::reactive({
    ai_model <- getUserOption(session, "llm_model")
    shiny::req(ai_model, ai_model != "")

    # Build context params for template substitution
    experiment <- pgx$name %||% pgx$description %||% "omics experiment"
    context_params <- list(
      experiment = experiment
    )

    omicsai::omicsai_config(
      model = ai_model,
      context = context_template,
      context_params = context_params
    )
  })

  # PlotModuleServer wrapper
  plot_server_wrapper <- function(id, func, func2, watermark) {
    PlotModuleServer(
      id,
      plotlib = "generic",
      plotlib2 = "generic",
      func = func,
      func2 = func2,
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI,
      pdf.width = 8, pdf.height = 5,
      res = c(75, 100),
      add.watermark = watermark
    )
  }

  # Initialize omicsai card module
  omicsai::omicsai_summary_card_server(
    id = id,
    params_reactive = params_reactive,
    template_reactive = template_reactive,
    config_reactive = config_reactive,
    plot_server_wrapper = plot_server_wrapper,
    cache = cache,
    watermark = watermark
  )
}
