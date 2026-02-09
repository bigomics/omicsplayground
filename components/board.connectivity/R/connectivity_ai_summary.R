##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for connectivity summary
#'
#' Extracts parameters from connectivity results to populate the prompt template.
#'
#' @param pgx PGX object
#' @param contrast Character; contrast name (e.g., "groupA_vs_groupB")
#' @param connectivity_scores Data frame; connectivity/similarity scores from
#'   getConnectivityScores() with columns: pathway, score, rho, NES, odd.ratio, tau
#' @param cum_enrichment_table Matrix; cumulative enrichment table from
#'   cumEnrichmentTable() (gene sets x signatures)
#' @param ntop Integer; number of top signatures to include (default 20)
#'
#' @return Named list with template parameters:
#'   contrast, top_signatures, leading_edge_genes, enriched_pathways, experiment
connectivity_build_ai_params <- function(pgx,
                                         contrast,
                                         connectivity_scores,
                                         cum_enrichment_table = NULL,
                                         ntop = 20) {
  fmt_num <- function(x, digits = 2) {
    ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
  }

  # Extract experiment description
  experiment <- ""
  if (!is.null(pgx$name) && nzchar(pgx$name)) {
    experiment <- pgx$name
  } else if (!is.null(pgx$description) && nzchar(pgx$description)) {
    experiment <- pgx$description
  }

  # Build top signatures section from connectivity_scores
  top_signatures <- ""
  if (!is.null(connectivity_scores) && is.data.frame(connectivity_scores) && nrow(connectivity_scores) > 0) {
    # Take top ntop by score (already sorted)
    tbl <- head(connectivity_scores, ntop)

    # Get available columns
    sig_cols <- intersect(c("pathway", "score", "rho", "NES", "odd.ratio", "tau"), colnames(tbl))

    sig_table <- paste0(
      "| Signature | Score | rho | NES | Odds ratio | tau |\n",
      "|-----------|-------|-----|-----|------------|-----|\n",
      paste(
        sprintf(
          "| %s | %s | %s | %s | %s | %s |",
          playbase::shortstring(tbl$pathway, 80),
          if ("score" %in% sig_cols) fmt_num(tbl$score, 3) else rep("NA", nrow(tbl)),
          if ("rho" %in% sig_cols) fmt_num(tbl$rho, 3) else rep("NA", nrow(tbl)),
          if ("NES" %in% sig_cols) fmt_num(tbl$NES, 3) else rep("NA", nrow(tbl)),
          if ("odd.ratio" %in% sig_cols) fmt_num(tbl$odd.ratio, 3) else rep("NA", nrow(tbl)),
          if ("tau" %in% sig_cols) fmt_num(tbl$tau, 3) else rep("NA", nrow(tbl))
        ),
        collapse = "\n"
      )
    )

    # Summary statistics
    n_concordant <- sum(tbl$rho > 0, na.rm = TRUE)
    n_discordant <- sum(tbl$rho < 0, na.rm = TRUE)

    top_signatures <- paste0(
      "**Top Similar Signatures (by connectivity score):**\n\n",
      sig_table, "\n\n",
      "**Direction:** ", n_concordant, " concordant (positive rho), ",
      n_discordant, " discordant (negative rho) of top ", nrow(tbl), " shown"
    )
  }

  # Build leading edge genes section
  leading_edge_genes <- ""
  if (!is.null(connectivity_scores) && is.data.frame(connectivity_scores) &&
    nrow(connectivity_scores) > 0 && "leadingEdge" %in% colnames(connectivity_scores)) {
    # Collect leading edge genes from top signatures
    top_sigs <- head(connectivity_scores, min(10, nrow(connectivity_scores)))
    le_genes <- unlist(top_sigs$leadingEdge)

    if (length(le_genes) > 0) {
      # Count frequency across signatures
      gene_freq <- sort(table(le_genes), decreasing = TRUE)
      top_le <- head(gene_freq, 15)

      # Get fold changes if available from pgx
      fc_vals <- rep(NA, length(top_le))
      if (!is.null(pgx$gx.meta$meta[[contrast]])) {
        mx <- pgx$gx.meta$meta[[contrast]]
        gene_probes <- NULL
        if (!is.null(pgx$genes) && "symbol" %in% colnames(pgx$genes)) {
          gene_probes <- rownames(pgx$genes)[pgx$genes$symbol %in% names(top_le)]
          gene_probes <- intersect(gene_probes, rownames(mx))
        }
        if (!is.null(gene_probes) && length(gene_probes) > 0) {
          probe_fc <- mx$meta.fx[gene_probes]
          probe_sym <- pgx$genes[gene_probes, "symbol"]
          fc_lookup <- setNames(probe_fc, probe_sym)
          fc_vals <- fc_lookup[names(top_le)]
        }
      }

      gene_table <- paste0(
        "| Gene | Frequency (across top signatures) | logFC |\n",
        "|------|-----------------------------------|-------|\n",
        paste(
          sprintf(
            "| %s | %s | %s |",
            names(top_le),
            as.integer(top_le),
            fmt_num(fc_vals, 3)
          ),
          collapse = "\n"
        )
      )

      leading_edge_genes <- paste0(
        "**Leading edge genes shared across top similar signatures (ranked by frequency):**\n\n",
        gene_table, "\n\n",
        "**Top ", length(top_le), " genes shown from the leading edges of the top ",
        nrow(top_sigs), " most similar signatures**"
      )
    }
  }

  if (leading_edge_genes == "") {
    leading_edge_genes <- "No leading edge gene information available for this analysis."
  }

  # Build enriched pathways section from cumulative enrichment table
  enriched_pathways <- ""
  if (!is.null(cum_enrichment_table) && is.matrix(cum_enrichment_table) &&
    nrow(cum_enrichment_table) > 0) {
    # Compute mean enrichment across signatures (excluding first column which is the query)
    if (ncol(cum_enrichment_table) > 1) {
      mean_enrich <- rowMeans(cum_enrichment_table[, -1, drop = FALSE], na.rm = TRUE)
    } else {
      mean_enrich <- cum_enrichment_table[, 1]
    }

    # Sort by absolute mean enrichment and take top
    top_idx <- head(order(-abs(mean_enrich)), min(15, length(mean_enrich)))
    top_pathways <- names(mean_enrich)[top_idx]
    top_values <- mean_enrich[top_idx]

    # Also get the query contrast enrichment
    query_values <- cum_enrichment_table[top_idx, 1]

    pathway_names <- sub(".*:", "", top_pathways)

    pw_table <- paste0(
      "| Pathway | Mean enrichment (across similar) | Query enrichment |\n",
      "|---------|----------------------------------|------------------|\n",
      paste(
        sprintf(
          "| %s | %s | %s |",
          pathway_names,
          fmt_num(top_values, 3),
          fmt_num(query_values, 3)
        ),
        collapse = "\n"
      )
    )

    n_up <- sum(top_values > 0, na.rm = TRUE)
    n_down <- sum(top_values < 0, na.rm = TRUE)

    enriched_pathways <- paste0(
      "**Top enriched pathways across similar signatures:**\n\n",
      pw_table, "\n\n",
      "**Direction:** ", n_up, " upregulated, ", n_down, " downregulated (of top ", length(top_idx), " shown)"
    )
  }

  if (enriched_pathways == "") {
    enriched_pathways <- "No cumulative enrichment data available for this analysis."
  }

  list(
    contrast = contrast,
    top_signatures = top_signatures,
    leading_edge_genes = leading_edge_genes,
    enriched_pathways = enriched_pathways,
    experiment = experiment
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: Connectivity AI Summary
# -----------------------------------------------------------------------------

#' Connectivity AI Summary UI
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
connectivity_ai_summary_ui <- function(id,
                                       title = "AI Summary",
                                       label = "",
                                       info.text = "",
                                       caption = "AI-generated connectivity summary.",
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

  omicsai::omicsai_summary_card_ui(
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

#' Connectivity AI Summary Server
#'
#' Server logic for connectivity AI-generated summaries.
#' Uses omicsai_summary_card_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param getConnectivityScores Reactive returning connectivity scores data frame
#' @param cumEnrichmentTable Reactive returning cumulative enrichment matrix
#' @param r_contrast Reactive returning selected contrast name
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_summary_card_server
connectivity_ai_summary_server <- function(id,
                                           pgx,
                                           getConnectivityScores,
                                           cumEnrichmentTable,
                                           r_contrast,
                                           session,
                                           cache = NULL,
                                           watermark = FALSE) {

  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    scores <- getConnectivityScores()
    contrast <- r_contrast()
    shiny::req(scores, contrast)

    # cumEnrichmentTable may be NULL
    cum_enrich <- tryCatch(cumEnrichmentTable(), error = function(e) NULL)

    connectivity_build_ai_params(
      pgx = pgx,
      contrast = contrast,
      connectivity_scores = scores,
      cum_enrichment_table = cum_enrich,
      ntop = 20
    )
  })

  # Load template from board.connectivity prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.connectivity/prompts/connectivity_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omicsai::omicsai_load_template(board_template_path)
  })

  # Load connectivity methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.connectivity/prompts/Connectivity_methods.md"
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
