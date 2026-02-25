##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for pathway enrichment summary
#'
#' Extracts top enriched pathways from the pathway table and formats
#' them as a markdown table for the AI prompt template.
#'
#' @param pgx PGX object containing experiment data
#' @param contrast Character; selected contrast name
#' @param pathway_table Data frame; pathway enrichment table with columns
#'   pathway, logFC, meta.q (from getReactomeTable, getWikiPathwayTable,
#'   or GO table_data)
#' @param pathway_type Character; pathway database type (e.g., "GO",
#'   "WikiPathways", "Reactome")
#' @param ntop Integer; number of top pathways to include (default 20)
#'
#' @return Named list with template parameters (board_params):
#'   contrast, pathway_type, genesets, experiment
pathway_build_ai_params <- function(pgx,
                                    contrast,
                                    pathway_table,
                                    pathway_type = "GO",
                                    ntop = 20) {

  # Build genesets section from pathway table
  genesets <- ""
  if (is.data.frame(pathway_table) && nrow(pathway_table) > 0) {
    # Determine column names based on pathway type
    # GO tables have: id, term, score, logFC, meta.q
    # Reactome/WikiPathway tables have: pathway (or pathway.id + pathway), logFC, meta.q
    has_term <- "term" %in% colnames(pathway_table)
    has_pathway <- "pathway" %in% colnames(pathway_table)
    has_score <- "score" %in% colnames(pathway_table)

    # Get pathway name column
    if (has_term) {
      pathway_names <- pathway_table$term
    } else if (has_pathway) {
      pathway_names <- pathway_table$pathway
    } else {
      pathway_names <- rownames(pathway_table)
    }

    # Clean pathway names (remove HTML links if present)
    pathway_names <- gsub("<[^>]+>", "", pathway_names)
    pathway_names <- trimws(pathway_names)

    # Get logFC values
    logfc <- if ("logFC" %in% colnames(pathway_table)) {
      pathway_table$logFC
    } else {
      rep(NA_real_, nrow(pathway_table))
    }

    # Get q-values
    qval <- if ("meta.q" %in% colnames(pathway_table)) {
      pathway_table$meta.q
    } else {
      rep(NA_real_, nrow(pathway_table))
    }

    # Sort by absolute logFC (most activated/repressed first)
    ord <- order(-abs(logfc))
    pathway_names <- pathway_names[ord]
    logfc <- logfc[ord]
    qval <- qval[ord]

    # Take top N
    n <- min(ntop, length(pathway_names))
    pathway_names <- pathway_names[seq_len(n)]
    logfc <- logfc[seq_len(n)]
    qval <- qval[seq_len(n)]

    # Build markdown table
    if (has_score) {
      score <- pathway_table$score[ord][seq_len(n)]
      gset_table <- paste0(
        "| Pathway | Score | logFC | Q-value |\n",
        "|---------|-------|-------|---------|\n",
        paste(
          sprintf(
            "| %s | %s | %s | %s |",
            pathway_names,
            omicsai::omicsai_format_num(score, 2),
            omicsai::omicsai_format_num(logfc, 2),
            omicsai::omicsai_format_num(qval, 4)
          ),
          collapse = "\n"
        )
      )
    } else {
      gset_table <- paste0(
        "| Pathway | logFC | Q-value |\n",
        "|---------|-------|---------|\n",
        paste(
          sprintf(
            "| %s | %s | %s |",
            pathway_names,
            omicsai::omicsai_format_num(logfc, 2),
            omicsai::omicsai_format_num(qval, 4)
          ),
          collapse = "\n"
        )
      )
    }

    # Summary statistics
    n_up <- sum(logfc > 0, na.rm = TRUE)
    n_down <- sum(logfc < 0, na.rm = TRUE)
    n_sig <- sum(qval < 0.05, na.rm = TRUE)

    genesets <- paste0(
      "**Top Enriched Pathways (by effect size):**\n\n",
      gset_table, "\n\n",
      "**Direction:** ", n_up, " upregulated, ", n_down, " downregulated\n",
      "**Significant (q < 0.05):** ", n_sig, " of ", n, " shown"
    )
  }

  # Get experiment description

  experiment <- pgx$description %||% pgx$name %||% ""

  list(
    contrast = contrast,
    pathway_type = pathway_type,
    genesets = genesets,
    experiment = experiment
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: Pathway AI Summary
# -----------------------------------------------------------------------------

#' Pathway AI Summary UI
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

#' Pathway AI Summary Server
#'
#' Server logic for pathway AI-generated enrichment summaries.
#' Uses omicsai_text_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (non-reactive) containing experiment data
#' @param pathway_table_reactive Reactive returning pathway enrichment data frame
#' @param gs_contrast Reactive returning selected contrast name
#' @param pathway_type Character; pathway database type (default "GO")
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_text_server
pathway_ai_summary_server <- function(id,
                                      pgx,
                                      pathway_table_reactive,
                                      gs_contrast,
                                      pathway_type = "GO",
                                      session,
                                      cache = NULL,
                                      watermark = FALSE) {

  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    tbl <- pathway_table_reactive()
    contrast <- gs_contrast()
    shiny::req(tbl, contrast)

    pathway_build_ai_params(
      pgx = pgx,
      contrast = contrast,
      pathway_table = tbl,
      pathway_type = pathway_type,
      ntop = 20
    )
  })

  # Load template from board.pathway prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.pathway/prompts/pathway_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omicsai::omicsai_load_template(board_template_path)
  })

  # Load Pathway methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.pathway/prompts/Pathway_methods.md"
  )
  context_template <- omicsai::omicsai_load_template(context_template_path)

  # Build config reactive from user's selected LLM model
  # Includes board-specific context in system prompt
  config_reactive <- shiny::reactive({
    ai_model <- getUserOption(session, "llm_model")
    shiny::req(ai_model, ai_model != "")

    # Build context params for template substitution
    context_params <- list(
      experiment = pgx$description %||% pgx$name %||% "omics experiment"
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
