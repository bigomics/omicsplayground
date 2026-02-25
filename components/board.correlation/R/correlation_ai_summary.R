##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for correlation summary
#'
#' Extracts parameters from correlation results to populate the prompt template.
#'
#' @param pgx PGX object
#' @param gene Character; selected gene (probe ID)
#' @param gene_corr Matrix; gene correlation matrix from getGeneCorr() with
#'   columns: cor, cov (and optionally others)
#' @param partial_corr Data frame; partial correlation from getPartialCorrelation()
#'   with columns: cor, pcor
#' @param gene_filter Character; gene filter applied (e.g., "<all>", gene family)
#' @param ntop Integer; number of top correlated genes to include (default 20)
#'
#' @return Named list with template parameters:
#'   gene, gene_symbol, gene_filter, correlated_genes, partial_correlation, experiment
correlation_build_ai_params <- function(pgx,
                                        gene,
                                        gene_corr,
                                        partial_corr,
                                        gene_filter = "<all>",
                                        ntop = 20) {

  # Extract experiment description
  experiment <- ""
  if (!is.null(pgx$name) && nzchar(pgx$name)) {
    experiment <- pgx$name
  } else if (!is.null(pgx$description) && nzchar(pgx$description)) {
    experiment <- pgx$description
  }

  # Get gene symbol
  gene_symbol <- gene
  if (!is.null(pgx$genes) && "symbol" %in% colnames(pgx$genes) &&
    gene %in% rownames(pgx$genes)) {
    sym <- pgx$genes[gene, "symbol"]
    if (!is.na(sym) && nzchar(sym)) gene_symbol <- sym
  }

  # Build correlated genes section from gene_corr
  correlated_genes <- ""
  if (!is.null(gene_corr) && is.matrix(gene_corr) && nrow(gene_corr) > 0) {
    # Remove the gene itself from the table
    tbl <- gene_corr[rownames(gene_corr) != gene, , drop = FALSE]
    # Sort by absolute correlation and take top ntop
    tbl <- tbl[order(-abs(tbl[, "cor"])), , drop = FALSE]
    tbl <- head(tbl, ntop)

    # Map probe IDs to symbols
    probes <- rownames(tbl)
    symbols <- probes
    if (!is.null(pgx$genes) && "symbol" %in% colnames(pgx$genes)) {
      matched <- intersect(probes, rownames(pgx$genes))
      if (length(matched) > 0) {
        sym_vals <- pgx$genes[matched, "symbol"]
        symbols[probes %in% matched] <- ifelse(
          is.na(sym_vals) | sym_vals == "", matched, sym_vals
        )
      }
    }

    # Get partial correlation values if available
    pcor_vals <- rep(NA, nrow(tbl))
    if (!is.null(partial_corr) && is.data.frame(partial_corr)) {
      matched_pcor <- intersect(probes, rownames(partial_corr))
      if (length(matched_pcor) > 0) {
        idx <- match(matched_pcor, probes)
        pcor_vals[idx] <- partial_corr[matched_pcor, "pcor"]
      }
    }

    # Build markdown table
    cor_table <- paste0(
      "| Gene | Correlation | Partial Corr | Covariance |\n",
      "|------|-------------|--------------|------------|\n",
      paste(
        sprintf(
          "| %s | %s | %s | %s |",
          symbols,
          omicsai::omicsai_format_num(tbl[, "cor"]),
          omicsai::omicsai_format_num(pcor_vals),
          omicsai::omicsai_format_num(tbl[, "cov"])
        ),
        collapse = "\n"
      )
    )

    # Summary statistics
    n_pos <- sum(tbl[, "cor"] > 0, na.rm = TRUE)
    n_neg <- sum(tbl[, "cor"] < 0, na.rm = TRUE)
    n_strong <- sum(abs(tbl[, "cor"]) > 0.6, na.rm = TRUE)

    correlated_genes <- paste0(
      "**Top Correlated Genes (by |correlation|):**\n\n",
      cor_table, "\n\n",
      "**Direction:** ", n_pos, " positively correlated, ", n_neg, " negatively correlated (of top ", nrow(tbl), " shown)  \n",
      "**Strongly correlated (|r| > 0.6):** ", n_strong, " of ", nrow(tbl), " shown"
    )
  }

  # Build partial correlation network section
  partial_correlation <- ""
  if (!is.null(partial_corr) && is.data.frame(partial_corr) && nrow(partial_corr) > 0) {
    # Remove the gene itself
    pcor_tbl <- partial_corr[rownames(partial_corr) != gene, , drop = FALSE]
    # Sort by absolute partial correlation
    pcor_tbl <- pcor_tbl[order(-abs(pcor_tbl$pcor)), , drop = FALSE]
    pcor_tbl <- head(pcor_tbl, 15)

    # Map probe IDs to symbols
    probes <- rownames(pcor_tbl)
    symbols <- probes
    if (!is.null(pgx$genes) && "symbol" %in% colnames(pgx$genes)) {
      matched <- intersect(probes, rownames(pgx$genes))
      if (length(matched) > 0) {
        sym_vals <- pgx$genes[matched, "symbol"]
        symbols[probes %in% matched] <- ifelse(
          is.na(sym_vals) | sym_vals == "", matched, sym_vals
        )
      }
    }

    pcor_table <- paste0(
      "| Gene | Partial Correlation | Correlation |\n",
      "|------|---------------------|-------------|\n",
      paste(
        sprintf(
          "| %s | %s | %s |",
          symbols,
          omicsai::omicsai_format_num(pcor_tbl$pcor),
          omicsai::omicsai_format_num(pcor_tbl$cor)
        ),
        collapse = "\n"
      )
    )

    n_direct_pos <- sum(pcor_tbl$pcor > 0.1, na.rm = TRUE)
    n_direct_neg <- sum(pcor_tbl$pcor < -0.1, na.rm = TRUE)

    partial_correlation <- paste0(
      "**Top Direct Interactions (by |partial correlation|):**\n\n",
      pcor_table, "\n\n",
      "**Direct associations (|pcor| > 0.1):** ", n_direct_pos, " positive, ", n_direct_neg, " negative (of top ", nrow(pcor_tbl), " shown)  \n",
      "**Top ", nrow(pcor_tbl), " genes shown ranked by absolute partial correlation**"
    )
  }

  list(
    gene = gene,
    gene_symbol = gene_symbol,
    gene_filter = gene_filter,
    correlated_genes = correlated_genes,
    partial_correlation = partial_correlation,
    experiment = experiment
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: Correlation AI Summary
# -----------------------------------------------------------------------------


#' Correlation AI Summary Server
#'
#' Server logic for correlation AI-generated summaries.
#' Uses omicsai_text_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param getGeneCorr Reactive returning gene correlation matrix
#' @param getPartialCorrelation Reactive returning partial correlation data frame
#' @param sel_gene Reactive returning selected gene (probe ID)
#' @param gene_filter Reactive returning gene filter selection
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_text_server
correlation_ai_summary_server <- function(id,
                                          pgx,
                                          getGeneCorr,
                                          getPartialCorrelation,
                                          sel_gene,
                                          gene_filter,
                                          session,
                                          cache = NULL,
                                          watermark = FALSE) {

  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    gene_corr <- getGeneCorr()
    partial_corr <- getPartialCorrelation()
    gene <- sel_gene()
    filter_val <- gene_filter()
    shiny::req(gene_corr, gene)

    correlation_build_ai_params(
      pgx = pgx,
      gene = gene,
      gene_corr = gene_corr,
      partial_corr = partial_corr,
      gene_filter = filter_val %||% "<all>",
      ntop = 20
    )
  })

  # Load template from board.correlation prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.correlation/prompts/correlation_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omicsai::omicsai_load_template(board_template_path)
  })

  # Load correlation methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.correlation/prompts/Correlation_methods.md"
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
      system_prompt = omicsai::omicsai_substitute_template(
        context_template, context_params
      )
    )
  })

  AiTextCardServer(
    id = id,
    params_reactive = params_reactive,
    template_reactive = template_reactive,
    config_reactive = config_reactive,
    cache = cache,
    watermark = watermark
  )
}
