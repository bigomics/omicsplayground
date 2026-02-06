##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for gene signature summary
#'
#' Extracts parameters from signature overlap and enrichment results
#' to populate the prompt template.
#'
#' @param pgx PGX object
#' @param markers List with features and symbols from getCurrentMarkers()
#' @param overlap_table Data frame from getOverlapTable() with columns:
#'   db, geneset, score, k/K, odds.ratio, q.fisher, common.genes
#' @param gene_table Data frame from getEnrichmentGeneTable() with columns:
#'   feature, symbol, log2FC, q.value (optional, may be NULL)
#' @param ntop Integer; number of top pathways/genes to include (default 20)
#'
#' @return Named list with template parameters:
#'   signature_size, genesets, top_genes, experiment
signature_build_ai_params <- function(pgx,
                                      markers,
                                      overlap_table,
                                      gene_table = NULL,
                                      ntop = 20) {
  fmt_num <- function(x, digits = 2) {
    ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
  }

  # Signature size
  signature_size <- length(markers$symbols)

  # Extract top pathways from overlap table
  genesets <- ""
  if (is.data.frame(overlap_table) && nrow(overlap_table) > 0) {
    # Already sorted by score in getOverlapTable()
    sig <- overlap_table[!is.na(overlap_table$q.fisher) &
      overlap_table$q.fisher < 0.05, , drop = FALSE]

    top_ov <- head(sig, ntop)
    if (nrow(top_ov) > 0) {
      pathway_names <- sub(".*:", "", top_ov$geneset)

      gset_table <- paste0(
        "| Pathway | Score | Odds Ratio | Q-value | Overlap (k/K) |\n",
        "|---------|-------|------------|---------|---------------|\n",
        paste(
          sprintf(
            "| %s | %s | %s | %s | %s |",
            pathway_names,
            fmt_num(top_ov$score, 2),
            fmt_num(top_ov$odds.ratio, 2),
            fmt_num(top_ov$q.fisher, 3),
            top_ov[["k/K"]]
          ),
          collapse = "\n"
        )
      )

      n_sig <- sum(overlap_table$q.fisher < 0.05, na.rm = TRUE)

      genesets <- paste0(
        "**Top Overlapping Gene Sets (q < 0.05):**\n\n",
        gset_table, "\n\n",
        "**Total significant overlaps:** ", n_sig, " of ",
        nrow(overlap_table), " tested"
      )
    }
  }

  if (genesets == "" && is.data.frame(overlap_table) && nrow(overlap_table) > 0) {
    # Fallback: just list top pathway names
    top_names <- head(sub(".*:", "", overlap_table$geneset), ntop)
    genesets <- paste0("- ", top_names, collapse = "\n")
  }

  # Extract top genes from gene_table if available
  top_genes <- ""
  if (is.data.frame(gene_table) && nrow(gene_table) > 0) {
    gene_table <- gene_table[order(-abs(gene_table$log2FC)), , drop = FALSE]
    top_g <- head(gene_table, ntop)

    gene_tbl <- paste0(
      "| Gene | Log2FC | Q-value |\n",
      "|------|--------|---------|\n",
      paste(
        sprintf(
          "| %s | %s | %s |",
          top_g$symbol,
          fmt_num(top_g$log2FC, 2),
          fmt_num(top_g$q.value, 3)
        ),
        collapse = "\n"
      )
    )

    n_up <- sum(top_g$log2FC > 0, na.rm = TRUE)
    n_down <- sum(top_g$log2FC < 0, na.rm = TRUE)

    top_genes <- paste0(
      "**Top Differentially Expressed Genes in Signature:**\n\n",
      gene_tbl, "\n\n",
      "**Top ", nrow(top_g), " of ", signature_size, " signature genes shown** ",
      "(", n_up, " up, ", n_down, " down)"
    )
  } else if (!is.null(markers$symbols) && length(markers$symbols) > 0) {
    # Fallback: just list gene symbols
    genes_list <- head(markers$symbols, ntop)
    top_genes <- paste0(
      "The following genes are in the signature:\n\n",
      paste(genes_list, collapse = ", ")
    )
  }

  # Get experiment description
  experiment <- pgx$description %||% pgx$name %||% ""

  list(
    signature_size = signature_size,
    genesets = genesets,
    top_genes = top_genes,
    experiment = experiment
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: Signature AI Summary
# -----------------------------------------------------------------------------

#' Signature AI Summary UI
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
signature_ai_summary_ui <- function(id,
                                    title = "AI Summary",
                                    label = "",
                                    info.text = "",
                                    caption = "AI-generated signature summary.",
                                    height = c("100%", TABLE_HEIGHT_MODAL),
                                    width = c("auto", "100%")) {
  # PlotModuleUI wrapper for omicsai card integration
  card_wrapper <- function(id, content, options, title, label, info.text,
                           caption, height, width, download.fmt, ...) {
    PlotModuleUI(
      id,
      outputFunc = shiny::htmlOutput,
      title = title,
      label = label,
      info.text = info.text,
      options = options,
      caption = caption,
      height = height,
      width = width,
      download.fmt = download.fmt
    )
  }

  omics.ai::omicsai_summary_card_ui(
    id = id,
    card_wrapper = card_wrapper,
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width
  )
}

#' Signature AI Summary Server
#'
#' Server logic for gene signature AI-generated summaries.
#' Uses omicsai_summary_card_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param getCurrentMarkers Reactive returning list(features, symbols)
#' @param getOverlapTable Reactive returning overlap data frame
#' @param getEnrichmentGeneTable Reactive returning gene-level data (optional)
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_summary_card_server
signature_ai_summary_server <- function(id,
                                        pgx,
                                        getCurrentMarkers,
                                        getOverlapTable,
                                        getEnrichmentGeneTable = shiny::reactive(NULL),
                                        session,
                                        cache = NULL,
                                        watermark = FALSE) {
  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    markers <- getCurrentMarkers()
    overlap <- getOverlapTable()
    gene_table <- getEnrichmentGeneTable()
    shiny::req(markers, overlap)

    signature_build_ai_params(
      pgx = pgx,
      markers = markers,
      overlap_table = overlap,
      gene_table = gene_table,
      ntop = 20
    )
  })

  # Load template from board.signature prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.signature/prompts/signature_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omics.ai::omicsai_load_template(board_template_path)
  })

  # Load Signature methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.signature/prompts/Signature_methods.md"
  )
  context_template <- omics.ai::omicsai_load_template(context_template_path)

  # Build config reactive from user's selected LLM model
  # Includes board-specific context in system prompt
  config_reactive <- shiny::reactive({
    ai_model <- getUserOption(session, "llm_model")
    shiny::req(ai_model, ai_model != "")

    # Build context params for template substitution
    context_params <- list(
      experiment = pgx$description %||% pgx$name %||% "omics experiment"
    )

    omics.ai::omicsai_config(
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
  omics.ai::omicsai_summary_card_server(
    id = id,
    params_reactive = params_reactive,
    template_reactive = template_reactive,
    config_reactive = config_reactive,
    plot_server_wrapper = plot_server_wrapper,
    cache = cache,
    watermark = watermark
  )
}
