##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for biomarker summary
#'
#' Extracts parameters from biomarker importance results to populate the prompt template.
#'
#' @param pgx PGX object
#' @param target Character; prediction target phenotype (column name from pgx$Y)
#' @param importance_result List; result from pgx.compute_importance() with $R matrix
#'   of importance scores (rows = features, cols = methods)
#' @param ntop Integer; number of top features to include (default 20)
#'
#' @return Named list with template parameters:
#'   target, phenotype, top_features, model_performance, experiment
biomarker_build_ai_params <- function(pgx,
                                      target,
                                      importance_result,
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

  # Extract phenotype from target name (readable form)
  phenotype <- gsub("_", " ", target)

  # Build top features section from importance result
  top_features <- ""
  if (!is.null(importance_result) && !is.null(importance_result$R)) {
    R <- importance_result$R

    # Sort by cumulative importance (sum across methods)
    cumulative <- rowSums(R, na.rm = TRUE)
    R <- R[order(-cumulative), , drop = FALSE]
    R <- head(R, ntop)

    # Convert probe IDs to gene symbols
    symbols <- playbase::probe2symbol(
      rownames(R), pgx$genes, "gene_name", fill_na = TRUE
    )

    # Build markdown table with per-method scores and cumulative rank
    methods <- colnames(R)
    header <- paste0(
      "| Gene | ", paste(methods, collapse = " | "), " | Cumulative |\n",
      "|------|", paste(rep("------|", length(methods)), collapse = ""), "------|\n"
    )

    rows <- vapply(seq_len(nrow(R)), function(i) {
      scores <- paste(fmt_num(R[i, ], 3), collapse = " | ")
      cum_score <- fmt_num(sum(R[i, ], na.rm = TRUE), 3)
      sprintf("| %s | %s | %s |", symbols[i], scores, cum_score)
    }, character(1))

    feature_table <- paste0(header, paste(rows, collapse = "\n"))

    # Summary statistics
    n_methods <- ncol(R)
    n_features <- nrow(R)

    top_features <- paste0(
      "**Top ", n_features, " Features by Cumulative Importance (across ", n_methods, " methods):**\n\n",
      feature_table, "\n\n",
      "**Methods used:** ", paste(methods, collapse = ", "), "  \n",
      "**Total features shown:** ", n_features, " of ", nrow(importance_result$R), " evaluated"
    )
  }

  # Build model performance section
  model_performance <- ""
  if (!is.null(importance_result)) {
    perf_parts <- c()

    # Extract decision tree info if available
    if (!is.null(importance_result$rf)) {
      tree_vars <- setdiff(importance_result$rf$frame$var, "<leaf>")
      if (length(tree_vars) > 0) {
        tree_symbols <- playbase::probe2symbol(
          tree_vars, pgx$genes, "gene_name", fill_na = TRUE
        )
        perf_parts <- c(perf_parts, paste0(
          "**Decision tree features:** ", paste(tree_symbols, collapse = ", ")
        ))
      }
    }

    # Number of classes in target
    if (!is.null(pgx$Y) && target %in% colnames(pgx$Y)) {
      classes <- unique(pgx$Y[, target])
      classes <- classes[!is.na(classes)]
      perf_parts <- c(perf_parts, paste0(
        "**Number of classes:** ", length(classes),
        " (", paste(head(classes, 10), collapse = ", "),
        if (length(classes) > 10) ", ..." else "", ")"
      ))
    }

    # Number of samples used
    if (!is.null(importance_result$X)) {
      perf_parts <- c(perf_parts, paste0(
        "**Samples used:** ", ncol(importance_result$X),
        "  \n**Features evaluated:** ", nrow(importance_result$X)
      ))
    }

    if (length(perf_parts) > 0) {
      model_performance <- paste(perf_parts, collapse = "  \n")
    } else {
      model_performance <- "Model performance metrics not available for this analysis."
    }
  }

  list(
    target = target,
    phenotype = phenotype,
    top_features = top_features,
    model_performance = model_performance,
    experiment = experiment
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: Biomarker AI Summary
# -----------------------------------------------------------------------------

#' Biomarker AI Summary UI
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
biomarker_ai_summary_ui <- function(id,
                                    title = "AI Summary",
                                    label = "",
                                    info.text = "",
                                    caption = "AI-generated biomarker summary.",
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

#' Biomarker AI Summary Server
#'
#' Server logic for biomarker AI-generated summaries.
#' Uses omicsai_summary_card_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param calcVariableImportance Reactive returning importance result from
#'   pgx.compute_importance() (with $R, $X, $rf components)
#' @param pdx_target Reactive returning selected prediction target
#' @param is_computed Reactive logical; TRUE after biomarker computation completes
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_summary_card_server
biomarker_ai_summary_server <- function(id,
                                        pgx,
                                        calcVariableImportance,
                                        pdx_target,
                                        is_computed,
                                        session,
                                        cache = NULL,
                                        watermark = FALSE) {
  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    shiny::req(is_computed())
    importance_result <- calcVariableImportance()
    target <- pdx_target()
    shiny::req(importance_result, target)

    biomarker_build_ai_params(
      pgx = pgx,
      target = target,
      importance_result = importance_result,
      ntop = 20
    )
  })

  # Load template from board.biomarker prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.biomarker/prompts/biomarker_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omics.ai::omicsai_load_template(board_template_path)
  })

  # Load biomarker methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.biomarker/prompts/Biomarker_methods.md"
  )
  context_template <- omics.ai::omicsai_load_template(context_template_path)

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
