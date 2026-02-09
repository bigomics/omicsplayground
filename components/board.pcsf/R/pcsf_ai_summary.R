##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for PCSF network summary
#'
#' Extracts parameters from PCSF network results to populate the prompt template.
#'
#' @param pgx PGX object
#' @param contrast Character; contrast name (e.g., "groupA_vs_groupB")
#' @param centrality_table Data frame; centrality table from pgx.getPCSFcentrality()
#'   with columns: symbol, gene_title, logFC, centrality
#' @param pcsf_graph igraph object; the PCSF solution graph
#' @param ntop Integer; number of top hub genes to include (default 15)
#'
#' @return Named list with template parameters:
#'   contrast, phenotype, network_summary, hub_genes, network_pathways, experiment
pcsf_build_ai_params <- function(pgx,
                                 contrast,
                                 centrality_table,
                                 pcsf_graph,
                                 ntop = 15) {

  # Extract experiment description
  experiment <- ""
  if (!is.null(pgx$name) && nzchar(pgx$name)) {
    experiment <- pgx$name
  } else if (!is.null(pgx$description) && nzchar(pgx$description)) {
    experiment <- pgx$description
  }

  # Extract phenotype from contrast name (readable form)
  phenotype <- gsub("_", " ", contrast)

  # Build network summary
  network_summary <- ""
  if (!is.null(pcsf_graph) && inherits(pcsf_graph, "igraph")) {
    n_nodes <- igraph::vcount(pcsf_graph)
    n_edges <- igraph::ecount(pcsf_graph)

    # Count terminal vs Steiner nodes
    node_types <- igraph::V(pcsf_graph)$type
    n_terminal <- sum(node_types == "Terminal", na.rm = TRUE)
    n_steiner <- sum(node_types == "Steiner", na.rm = TRUE)

    # Get fold changes for network nodes
    F <- playbase::pgx.getMetaMatrix(pgx, level = "gene")$fc
    node_names <- igraph::V(pcsf_graph)$name
    fc_vals <- F[intersect(node_names, rownames(F)), contrast]
    fc_vals <- fc_vals[!is.na(fc_vals)]

    n_up <- sum(fc_vals > 0, na.rm = TRUE)
    n_down <- sum(fc_vals < 0, na.rm = TRUE)

    network_summary <- paste0(
      "**Network Topology:**\n",
      "- **Total nodes:** ", n_nodes, "\n",
      "- **Total edges:** ", n_edges, "\n",
      "- **Terminal nodes (DE genes):** ", n_terminal, "\n",
      "- **Steiner nodes (connectors):** ", n_steiner, "\n",
      "- **Upregulated nodes:** ", n_up, "\n",
      "- **Downregulated nodes:** ", n_down
    )
  }

  # Build hub genes section from centrality table
  hub_genes <- ""
  if (!is.null(centrality_table) && is.data.frame(centrality_table) && nrow(centrality_table) > 0) {
    tbl <- centrality_table
    tbl <- tbl[order(-tbl$centrality), , drop = FALSE]
    tbl <- head(tbl, ntop)

    # Get symbol column
    sym_col <- if ("symbol" %in% colnames(tbl)) tbl$symbol else rownames(tbl)

    # Get gene_title if available
    title_col <- if ("gene_title" %in% colnames(tbl)) tbl$gene_title else rep("", nrow(tbl))
    title_col[is.na(title_col)] <- ""

    hub_table <- paste0(
      "| Gene | logFC | Centrality | Description |\n",
      "|------|-------|------------|-------------|\n",
      paste(
        sprintf(
          "| %s | %s | %s | %s |",
          sym_col,
          omicsai::format_num(tbl$logFC, 3),
          omicsai::format_num(tbl$centrality, 4),
          substr(title_col, 1, 50)
        ),
        collapse = "\n"
      )
    )

    hub_genes <- paste0(
      "**Top Hub Genes (ranked by centrality score):**\n\n",
      hub_table, "\n\n",
      "**Top ", nrow(tbl), " genes shown by network centrality**"
    )
  }

  # Build network pathways section
  network_pathways <- ""
  if (!is.null(pcsf_graph) && inherits(pcsf_graph, "igraph")) {
    node_names <- igraph::V(pcsf_graph)$name

    # Try to get geneset-level enrichment for network genes
    if (!is.null(pgx$gsetX) && !is.null(pgx$GMT)) {
      # Find gene sets that overlap with network nodes
      gmt_genes <- rownames(pgx$GMT)
      network_genes <- intersect(node_names, gmt_genes)

      if (length(network_genes) > 5) {
        # Get geneset fold changes for this contrast
        gs_fc <- playbase::pgx.getMetaMatrix(pgx, level = "geneset")$fc
        if (!is.null(gs_fc) && contrast %in% colnames(gs_fc)) {
          gs_vals <- gs_fc[, contrast]
          gs_vals <- gs_vals[!is.na(gs_vals)]
          gs_vals <- sort(gs_vals, decreasing = TRUE)

          # Filter to gene sets with network overlap
          gmt_matrix <- pgx$GMT[network_genes, , drop = FALSE]
          overlap_counts <- colSums(gmt_matrix != 0)
          relevant_gsets <- names(overlap_counts[overlap_counts >= 3])
          relevant_gsets <- intersect(relevant_gsets, names(gs_vals))

          if (length(relevant_gsets) > 0) {
            gs_sub <- gs_vals[relevant_gsets]
            gs_sub <- gs_sub[order(-abs(gs_sub))]
            gs_sub <- head(gs_sub, 10)

            pathway_names <- sub(".*:", "", names(gs_sub))
            overlap_n <- overlap_counts[names(gs_sub)]

            pathway_table <- paste0(
              "| Pathway | logFC | Network Genes |\n",
              "|---------|-------|---------------|\n",
              paste(
                sprintf(
                  "| %s | %s | %s |",
                  pathway_names,
                  omicsai::format_num(gs_sub, 3),
                  overlap_n
                ),
                collapse = "\n"
              )
            )

            network_pathways <- paste0(
              "**Pathways enriched among PCSF network genes:**\n\n",
              pathway_table, "\n\n",
              "**Top ", length(gs_sub), " pathways shown (with >= 3 network genes)**"
            )
          }
        }
      }
    }
  }

  if (network_pathways == "") {
    network_pathways <- "No pathway enrichment data available for network genes."
  }

  list(
    contrast = contrast,
    phenotype = phenotype,
    network_summary = network_summary,
    hub_genes = hub_genes,
    network_pathways = network_pathways,
    experiment = experiment
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: PCSF AI Summary
# -----------------------------------------------------------------------------

#' PCSF AI Summary UI
#'
#' Creates AI summary UI wrapped in PlotModuleUI for consistent OmicsPlayground styling.
#'
#' @param id Shiny module namespace ID
#' @param title Card title
#' @param label Optional label
#' @param info.text Help/info text
#' @param caption Caption text
#' @param height Card height (can be vector for responsive sizing)
#' @param width Card width (can be vector for responsive sizing)
#'
#' @return Shiny UI element

#' PCSF AI Summary Server
#'
#' Server logic for PCSF AI-generated network summaries.
#' Uses omicsai_summary_card_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param getCentralityTable Reactive returning centrality data frame
#' @param getPcsfGraph Reactive returning PCSF igraph object
#' @param r_contrast Reactive returning selected contrast
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_summary_card_server
pcsf_ai_summary_server <- function(id,
                                   pgx,
                                   getCentralityTable,
                                   getPcsfGraph,
                                   r_contrast,
                                   session,
                                   cache = NULL,
                                   watermark = FALSE) {

  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    centrality_table <- getCentralityTable()
    pcsf_graph <- getPcsfGraph()
    contrast <- r_contrast()
    shiny::req(centrality_table, pcsf_graph, contrast)

    pcsf_build_ai_params(
      pgx = pgx,
      contrast = contrast,
      centrality_table = centrality_table,
      pcsf_graph = pcsf_graph,
      ntop = 15
    )
  })

  # Load template from board.pcsf prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.pcsf/prompts/pcsf_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omicsai::omicsai_load_template(board_template_path)
  })

  # Load PCSF methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.pcsf/prompts/PCSF_methods.md"
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
