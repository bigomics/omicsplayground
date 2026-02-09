##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for drug connectivity summary
#'
#' Extracts parameters from drug connectivity / DSEA results to populate
#' the prompt template.
#'
#' @param pgx PGX object
#' @param contrast Character; contrast name (e.g., "groupA_vs_groupB")
#' @param method Character; drug enrichment method name
#' @param dsea_table Data frame; drug enrichment table from getActiveDSEA()
#'   with columns: drug, NES, pval, padj, moa, target
#' @param moa_class Data frame; MOA class enrichment results from getMOA.class()
#' @param moa_target Data frame; MOA target enrichment results from getMOA.target()
#' @param ntop Integer; number of top drugs to include (default 20)
#'
#' @return Named list with template parameters:
#'   contrast, method, top_drugs, moa_summary, experiment
drugconnectivity_build_ai_params <- function(pgx,
                                             contrast,
                                             method,
                                             dsea_table,
                                             moa_class = NULL,
                                             moa_target = NULL,
                                             ntop = 20) {

  # Extract experiment description
  experiment <- pgx$name %||% pgx$description %||% "omics experiment"

  # Build top drugs section from dsea_table
  top_drugs <- ""
  if (!is.null(dsea_table) && is.data.frame(dsea_table) && nrow(dsea_table) > 0) {
    # Take top ntop by absolute NES
    tbl <- dsea_table
    tbl <- tbl[order(-abs(tbl$NES)), , drop = FALSE]
    tbl <- head(tbl, ntop)

    drug_table <- paste0(
      "| Drug | NES | p-value | q-value | MOA | Target |\n",
      "|------|-----|---------|---------|-----|--------|\n",
      paste(
        sprintf(
          "| %s | %s | %s | %s | %s | %s |",
          tbl$drug,
          omicsai::format_num(tbl$NES, 3),
          omicsai::format_num(tbl$pval, 4),
          omicsai::format_num(tbl$padj, 4),
          ifelse(is.na(tbl$moa) | tbl$moa == "", "-", tbl$moa),
          ifelse(is.na(tbl$target) | tbl$target == "", "-", tbl$target)
        ),
        collapse = "\n"
      )
    )

    # Summary statistics
    n_positive <- sum(tbl$NES > 0, na.rm = TRUE)
    n_negative <- sum(tbl$NES < 0, na.rm = TRUE)
    n_sig <- sum(dsea_table$padj < 0.05, na.rm = TRUE)

    top_drugs <- paste0(
      "**Top Connected Drugs (by |NES|):**\n\n",
      drug_table, "\n\n",
      "**Direction:** ", n_positive, " positive (mimicking), ",
      n_negative, " negative (opposing) connectivity (of top ", nrow(tbl), " shown)  \n",
      "**Total significant drugs (q < 0.05):** ", n_sig, " of ", nrow(dsea_table), " tested"
    )
  }

  # Build MOA summary section
  moa_summary <- ""
  moa_parts <- c()

  # MOA class section
  if (!is.null(moa_class) && is.data.frame(moa_class) && nrow(moa_class) > 0) {
    mc <- head(moa_class[order(-abs(moa_class$NES)), , drop = FALSE], 15)
    moa_class_table <- paste0(
      "**Top MOA classes:**\n\n",
      "| MOA Class | NES | p-value | q-value | Size |\n",
      "|-----------|-----|---------|---------|------|\n",
      paste(
        sprintf(
          "| %s | %s | %s | %s | %s |",
          mc$pathway,
          omicsai::format_num(mc$NES, 3),
          omicsai::format_num(mc$pval, 4),
          omicsai::format_num(mc$padj, 4),
          mc$size
        ),
        collapse = "\n"
      )
    )
    moa_parts <- c(moa_parts, moa_class_table)
  }

  # MOA target section
  if (!is.null(moa_target) && is.data.frame(moa_target) && nrow(moa_target) > 0) {
    mt <- head(moa_target[order(-abs(moa_target$NES)), , drop = FALSE], 15)
    moa_target_table <- paste0(
      "**Top molecular targets:**\n\n",
      "| Target | NES | p-value | q-value | Size |\n",
      "|--------|-----|---------|---------|------|\n",
      paste(
        sprintf(
          "| %s | %s | %s | %s | %s |",
          mt$pathway,
          omicsai::format_num(mt$NES, 3),
          omicsai::format_num(mt$pval, 4),
          omicsai::format_num(mt$padj, 4),
          mt$size
        ),
        collapse = "\n"
      )
    )
    moa_parts <- c(moa_parts, moa_target_table)
  }

  if (length(moa_parts) > 0) {
    moa_summary <- paste(moa_parts, collapse = "\n\n")
  }

  list(
    contrast = contrast,
    method = method,
    top_drugs = top_drugs,
    moa_summary = moa_summary,
    experiment = experiment
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: Drug Connectivity AI Summary
# -----------------------------------------------------------------------------

#' Drug Connectivity AI Summary UI
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

#' Drug Connectivity AI Summary Server
#'
#' Server logic for drug connectivity AI-generated summaries.
#' Uses omicsai_summary_card_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param getActiveDSEA Reactive returning active DSEA results (list with $table)
#' @param getMOA.class Reactive returning MOA class enrichment results
#' @param getMOA.target Reactive returning MOA target enrichment results
#' @param dsea_contrast Reactive returning selected contrast
#' @param dsea_method Reactive returning selected analysis method
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_summary_card_server
drugconnectivity_ai_summary_server <- function(id,
                                               pgx,
                                               getActiveDSEA,
                                               getMOA.class,
                                               getMOA.target,
                                               dsea_contrast,
                                               dsea_method,
                                               session,
                                               cache = NULL,
                                               watermark = FALSE) {

  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    dsea <- getActiveDSEA()
    contrast <- dsea_contrast()
    method <- dsea_method()
    shiny::req(dsea, contrast, method)

    dsea_table <- dsea$table

    # Safely attempt to get MOA data (may fail if no annotated drugs)
    moa_class <- tryCatch(getMOA.class(), error = function(e) NULL)
    moa_target <- tryCatch(getMOA.target(), error = function(e) NULL)

    drugconnectivity_build_ai_params(
      pgx = pgx,
      contrast = contrast,
      method = method,
      dsea_table = dsea_table,
      moa_class = moa_class,
      moa_target = moa_target,
      ntop = 20
    )
  })

  # Load template from board.drugconnectivity prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.drugconnectivity/prompts/drugconnectivity_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omicsai::omicsai_load_template(board_template_path)
  })

  # Load drug connectivity methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.drugconnectivity/prompts/DrugConnectivity_methods.md"
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
