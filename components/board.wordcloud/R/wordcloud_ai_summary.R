##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for WordCloud keyword enrichment summary
#'
#' Extracts parameters from word enrichment results to populate the prompt template.
#'
#' @param word_enrichment Data frame from getCurrentWordEnrichment();
#'   columns: word, NES, pval, padj, size, ES
#' @param contrast Character; current contrast name
#' @param ntop Integer; number of top keywords to include (default 30)
#'
#' @return Named list with template parameters:
#'   contrast, top_keywords, keyword_stats, experiment
wordcloud_build_ai_params <- function(word_enrichment,
                                      contrast,
                                      ntop = 30) {

  # Select top keywords by absolute NES
  df <- word_enrichment
  df <- df[order(-abs(df$NES)), , drop = FALSE]
  top_df <- head(df, ntop)

  # Format as markdown table
  top_keywords <- paste0(
    "| Keyword | NES | padj | Size |\n",
    "|---------|-----|------|------|\n",
    paste(
      sprintf(
        "| %s | %s | %s | %s |",
        top_df$word,
        omicsai::format_num(top_df$NES, 2),
        omicsai::format_num(top_df$padj, 3),
        top_df$size
      ),
      collapse = "\n"
    )
  )

  # Compute summary statistics
  n_sig <- sum(df$padj < 0.05, na.rm = TRUE)
  n_total <- nrow(df)
  n_pos <- sum(df$NES > 0 & df$padj < 0.05, na.rm = TRUE)
  n_neg <- sum(df$NES < 0 & df$padj < 0.05, na.rm = TRUE)

  nes_vals <- df$NES[!is.na(df$NES)]
  nes_range <- if (length(nes_vals) > 0) {
    paste0(
      omicsai::format_num(min(nes_vals), 2), " to ",
      omicsai::format_num(max(nes_vals), 2),
      " (median: ", omicsai::format_num(median(nes_vals), 2), ")"
    )
  } else {
    "NA"
  }

  keyword_stats <- paste0(
    "**Keyword Enrichment Statistics:**\n",
    "- **Significant keywords (padj < 0.05):** ", n_sig, " of ", n_total, " tested\n",
    "- **Positively enriched:** ", n_pos, " keywords\n",
    "- **Negatively enriched:** ", n_neg, " keywords\n",
    "- **NES range:** ", nes_range
  )

  list(
    contrast = contrast,
    top_keywords = top_keywords,
    keyword_stats = keyword_stats,
    experiment = ""
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: WordCloud AI Summary
# -----------------------------------------------------------------------------

#' WordCloud AI Summary UI
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

#' WordCloud AI Summary Server
#'
#' Server logic for WordCloud AI-generated keyword enrichment summaries.
#' Uses omicsai_summary_card_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (used for experiment description)
#' @param getCurrentWordEnrichment Reactive returning word enrichment data frame
#' @param r_contrast Reactive returning current contrast name
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_summary_card_server
wordcloud_ai_summary_server <- function(id,
                                        pgx,
                                        getCurrentWordEnrichment,
                                        r_contrast,
                                        session,
                                        cache = NULL,
                                        watermark = FALSE) {

  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    we <- getCurrentWordEnrichment()
    contrast <- r_contrast()
    shiny::req(we, contrast)

    wordcloud_build_ai_params(
      word_enrichment = we,
      contrast = contrast,
      ntop = 30
    )
  })

  # Load template from board.wordcloud prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.wordcloud/prompts/wordcloud_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omicsai::omicsai_load_template(board_template_path)
  })

  # Load WordCloud methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.wordcloud/prompts/WordCloud_methods.md"
  )
  context_template <- omicsai::omicsai_load_template(context_template_path)

  # Build config reactive from user's selected LLM model
  # Includes board-specific context in system prompt
  config_reactive <- shiny::reactive({
    ai_model <- getUserOption(session, "llm_model")
    shiny::req(ai_model, ai_model != "")

    # Build context params for template substitution
    context_params <- list(
      experiment = pgx$experiment %||% "omics experiment"
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
