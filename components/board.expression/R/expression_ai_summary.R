##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for expression summary
#'
#' Extracts parameters from differential expression results to populate the
#' prompt template.
#'
#' @param pgx PGX object
#' @param contrast Character; contrast name (e.g., "groupA_vs_groupB")
#' @param ntop Integer; number of top DE genes to include (default 20)
#'
#' @return Named list with template parameters:
#'   contrast, phenotype, top_genes, summary_stats, experiment
expression_build_ai_params <- function(pgx, contrast, ntop = 20) {
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

  # Extract phenotype from contrast name (readable form)
  phenotype <- gsub("_", " ", contrast)

  # Get gene-level meta results for this contrast
  mx <- pgx$gx.meta$meta[[contrast]]

  top_genes <- ""
  summary_stats <- ""

  if (!is.null(mx)) {
    # Compute meta fold change and meta q-value
    meta_fx <- mx$meta.fx
    meta_q <- mx$meta.q
    names(meta_fx) <- rownames(mx)
    names(meta_q) <- rownames(mx)

    # Remove NA / Inf values
    valid <- !is.na(meta_fx) & !is.infinite(meta_fx) & !is.na(meta_q)
    meta_fx <- meta_fx[valid]
    meta_q <- meta_q[valid]

    # Get gene symbols
    symbols <- rep(NA, length(meta_fx))
    if (!is.null(pgx$genes) && "symbol" %in% colnames(pgx$genes)) {
      symbols <- pgx$genes[names(meta_fx), "symbol"]
      symbols <- ifelse(is.na(symbols) | symbols == "", names(meta_fx), symbols)
    } else {
      symbols <- names(meta_fx)
    }

    # Build data frame and sort by absolute fold change
    df <- data.frame(
      probe = names(meta_fx),
      symbol = symbols,
      logFC = meta_fx,
      meta.q = meta_q,
      stringsAsFactors = FALSE
    )
    df <- df[order(-abs(df$logFC)), , drop = FALSE]

    # Deduplicate by symbol, keeping highest absolute FC
    df <- df[!duplicated(df$symbol), , drop = FALSE]

    # Summary statistics (before subsetting to ntop)
    n_total <- nrow(df)
    n_up <- sum(df$logFC > 0 & df$meta.q < 0.05, na.rm = TRUE)
    n_down <- sum(df$logFC < 0 & df$meta.q < 0.05, na.rm = TRUE)
    n_sig <- n_up + n_down

    summary_stats <- paste0(
      "- **Total genes tested:** ", n_total, "\n",
      "- **Significant genes (q < 0.05):** ", n_sig, "\n",
      "  - Upregulated: ", n_up, "\n",
      "  - Downregulated: ", n_down, "\n",
      "- **Top ", min(ntop, nrow(df)), " genes shown** (ranked by |logFC|)"
    )

    # Take top ntop genes
    df <- head(df, ntop)

    if (nrow(df) > 0) {
      gene_table <- paste0(
        "| Gene | logFC | q-value |\n",
        "|------|-------|---------|\n",
        paste(
          sprintf(
            "| %s | %s | %s |",
            df$symbol,
            fmt_num(df$logFC, 3),
            fmt_num(df$meta.q, 4)
          ),
          collapse = "\n"
        )
      )

      n_up_top <- sum(df$logFC > 0, na.rm = TRUE)
      n_down_top <- sum(df$logFC < 0, na.rm = TRUE)

      top_genes <- paste0(
        "**Top Differentially Expressed Genes (by |logFC|):**\n\n",
        gene_table, "\n\n",
        "**Direction:** ", n_up_top, " upregulated, ", n_down_top,
        " downregulated (of top ", nrow(df), " shown)"
      )
    }
  }

  list(
    contrast = contrast,
    phenotype = phenotype,
    top_genes = top_genes,
    summary_stats = summary_stats,
    experiment = experiment
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: Expression AI Summary
# -----------------------------------------------------------------------------

#' Expression AI Summary UI
#'
#' Creates AI summary UI wrapped in PlotModuleUI for consistent OmicsPlayground
#' styling.
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
expression_ai_summary_ui <- function(id,
                                     title = "AI Summary",
                                     label = "",
                                     info.text = "",
                                     caption = "AI-generated differential expression summary.",
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

#' Expression AI Summary Server
#'
#' Server logic for expression AI-generated summaries.
#' Uses omicsai_summary_card_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param contrast_reactive Reactive returning selected contrast
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_summary_card_server
expression_ai_summary_server <- function(id,
                                         pgx,
                                         contrast_reactive,
                                         session,
                                         cache = NULL,
                                         watermark = FALSE) {

  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    contrast <- contrast_reactive()
    shiny::req(contrast)

    expression_build_ai_params(
      pgx = pgx,
      contrast = contrast,
      ntop = 20
    )
  })

  # Load template from board.expression prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.expression/prompts/expression_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omicsai::omicsai_load_template(board_template_path)
  })

  # Load expression methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.expression/prompts/Expression_methods.md"
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
