##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for enrichment summary
#'
#' Extracts parameters from enrichment results to populate the prompt template.
#'
#' @param pgx PGX object
#' @param contrast Character; contrast name (e.g., "groupA_vs_groupB")
#' @param filtered_table Data frame; filtered gene set enrichment table from
#'   getFilteredGeneSetTable() with columns: logFC, meta.q, matched.genes/size, stars
#' @param ntop Integer; number of top gene sets to include (default 20)
#'
#' @return Named list with template parameters:
#'   contrast, phenotype, genesets, top_genes, experiment
enrichment_build_ai_params <- function(pgx,
                                       contrast,
                                       filtered_table,
                                       ntop = 20) {

  # Extract experiment description
  experiment <- ""
  if (!is.null(pgx$name) && nzchar(pgx$name)) {
    experiment <- pgx$name
  } else if (!is.null(pgx$description) && nzchar(pgx$description)) {
    experiment <- pgx$description
  }

  # Extract phenotype from contrast name (readable form)
  phenotype <- gsub("_", " ", contrast)

  # Build gene sets section from filtered_table
  genesets <- ""
  if (!is.null(filtered_table) && is.data.frame(filtered_table) && nrow(filtered_table) > 0) {
    # Sort by absolute logFC and take top ntop
    tbl <- filtered_table
    tbl <- tbl[order(-abs(tbl$logFC)), , drop = FALSE]
    tbl <- head(tbl, ntop)

    # Extract pathway names (strip collection prefix if present)
    pathways <- sub(".*:", "", rownames(tbl))

    # Get size column (matched.genes or size)
    size_col <- intersect(c("matched.genes", "size"), colnames(tbl))
    sizes <- if (length(size_col) > 0) tbl[[size_col[1]]] else rep(NA, nrow(tbl))

    # Get stars column
    stars_col <- if ("stars" %in% colnames(tbl)) tbl$stars else rep("", nrow(tbl))

    gset_table <- paste0(
      "| Pathway | logFC | q-value | Size | Stars |\n",
      "|---------|-------|---------|------|-------|\n",
      paste(
        sprintf(
          "| %s | %s | %s | %s | %s |",
          pathways,
          omicsai::format_num(tbl$logFC, 3),
          omicsai::format_num(tbl$meta.q, 4),
          sizes,
          stars_col
        ),
        collapse = "\n"
      )
    )

    # Summary statistics
    n_up <- sum(tbl$logFC > 0, na.rm = TRUE)
    n_down <- sum(tbl$logFC < 0, na.rm = TRUE)
    n_sig <- sum(filtered_table$meta.q < 0.05, na.rm = TRUE)

    genesets <- paste0(
      "**Top Enriched Gene Sets (by |logFC|):**\n\n",
      gset_table, "\n\n",
      "**Direction:** ", n_up, " upregulated, ", n_down, " downregulated (of top ", nrow(tbl), " shown)  \n",
      "**Total significant gene sets (q < 0.05):** ", n_sig, " of ", nrow(filtered_table), " tested"
    )
  }

  # Build leading edge / top genes section
  top_genes <- ""
  if (!is.null(filtered_table) && is.data.frame(filtered_table) && nrow(filtered_table) > 0) {
    # Get gene-level info from pgx for the top gene sets
    top_gsets <- head(rownames(filtered_table[order(-abs(filtered_table$logFC)), , drop = FALSE]), 5)

    gene_list <- c()
    for (gs in top_gsets) {
      if (!is.null(pgx$GMT) && gs %in% colnames(pgx$GMT)) {
        genes_in_set <- rownames(pgx$GMT)[which(pgx$GMT[, gs] != 0)]
        # Map to symbols if possible
        if (!is.null(pgx$genes) && "symbol" %in% colnames(pgx$genes)) {
          matched <- intersect(genes_in_set, rownames(pgx$genes))
          if (length(matched) > 0) {
            symbols <- pgx$genes[matched, "symbol"]
            symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
            gene_list <- c(gene_list, symbols)
          }
        } else {
          gene_list <- c(gene_list, genes_in_set)
        }
      }
    }

    # Get fold changes for these genes from the contrast
    if (length(gene_list) > 0 && !is.null(pgx$gx.meta$meta[[contrast]])) {
      mx <- pgx$gx.meta$meta[[contrast]]
      gene_probes <- NULL
      if (!is.null(pgx$genes) && "symbol" %in% colnames(pgx$genes)) {
        gene_probes <- rownames(pgx$genes)[pgx$genes$symbol %in% gene_list]
        gene_probes <- intersect(gene_probes, rownames(mx))
      }

      if (!is.null(gene_probes) && length(gene_probes) > 0) {
        fc_vals <- mx$meta.fx[gene_probes]
        symbols <- pgx$genes[gene_probes, "symbol"]
        symbols <- ifelse(is.na(symbols) | symbols == "", gene_probes, symbols)

        # Deduplicate by symbol, keeping highest abs FC
        df <- data.frame(
          symbol = symbols,
          fc = fc_vals,
          stringsAsFactors = FALSE
        )
        df <- df[!is.na(df$fc), , drop = FALSE]
        df <- df[order(-abs(df$fc)), , drop = FALSE]
        df <- df[!duplicated(df$symbol), , drop = FALSE]
        df <- head(df, 15)

        if (nrow(df) > 0) {
          gene_table <- paste0(
            "| Gene | logFC |\n",
            "|------|-------|\n",
            paste(
              sprintf("| %s | %s |", df$symbol, omicsai::format_num(df$fc, 3)),
              collapse = "\n"
            )
          )

          top_genes <- paste0(
            "**Top genes from leading enriched gene sets (ranked by |logFC|):**\n\n",
            gene_table, "\n\n",
            "**Top ", nrow(df), " genes shown from the top ", length(top_gsets), " enriched gene sets**"
          )
        }
      }
    }

    if (top_genes == "") {
      # Fallback: just list unique gene names
      gene_list <- unique(gene_list)
      if (length(gene_list) > 0) {
        gene_list <- head(gene_list, 30)
        top_genes <- paste0(
          "Genes frequently appearing in the top enriched gene sets:\n\n",
          paste(gene_list, collapse = ", ")
        )
      }
    }
  }

  list(
    contrast = contrast,
    phenotype = phenotype,
    genesets = genesets,
    top_genes = top_genes,
    experiment = experiment
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: Enrichment AI Summary
# -----------------------------------------------------------------------------

#' Enrichment AI Summary UI
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

#' Enrichment AI Summary Server
#'
#' Server logic for enrichment AI-generated summaries.
#' Uses omicsai_summary_card_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param getFilteredGeneSetTable Reactive returning filtered gene set table
#' @param gs_contrast Reactive returning selected contrast
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_summary_card_server
enrichment_ai_summary_server <- function(id,
                                         pgx,
                                         getFilteredGeneSetTable,
                                         gs_contrast,
                                         session,
                                         cache = NULL,
                                         watermark = FALSE) {

  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    filtered_table <- getFilteredGeneSetTable()
    contrast <- gs_contrast()
    shiny::req(filtered_table, contrast)

    enrichment_build_ai_params(
      pgx = pgx,
      contrast = contrast,
      filtered_table = filtered_table,
      ntop = 20
    )
  })

  # Load template from board.enrichment prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.enrichment/prompts/enrichment_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omicsai::omicsai_load_template(board_template_path)
  })

  # Load enrichment methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.enrichment/prompts/Enrichment_methods.md"
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
