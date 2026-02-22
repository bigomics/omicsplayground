##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for time series cluster summary
#'
#' Extracts parameters from time series clustering results to populate the
#' prompt template. Formats enrichment data (gene set correlations) and
#' cluster statistics for the selected module.
#'
#' @param data List from timeseries_filtered() reactive with Z, X, modules, gset.rho
#' @param pgx PGX object with gene annotations and metadata
#' @param module Character; selected module name (e.g., "M1", "M2")
#' @param contrast Character; selected contrast name
#' @param ntop Integer; number of top gene sets to include (default 20)
#'
#' @return Named list with template parameters:
#'   contrast, genesets, cluster_info, experiment
timeseries_build_ai_params <- function(data,
                                       pgx,
                                       module,
                                       contrast,
                                       ntop = 20) {

  ## -----------------------------------------------------------
  ## Cluster information
  ## -----------------------------------------------------------
  module_genes <- names(data$modules)[data$modules == module]
  n_genes <- length(module_genes)

  ## Compute module temporal profile summary
  cluster_lines <- c("**Cluster Statistics:**")
  cluster_lines <- c(cluster_lines, paste0("- **Size:** ", n_genes, " genes"))

  if (!is.null(data$Z) && length(module_genes) > 0) {
    module_expr <- data$Z[module_genes, , drop = FALSE]
    avg_profile <- colMeans(module_expr, na.rm = TRUE)
    time_points <- colnames(data$Z)

    ## Determine temporal trend direction
    if (length(avg_profile) >= 2) {
      first_half <- mean(avg_profile[1:ceiling(length(avg_profile) / 2)], na.rm = TRUE)
      second_half <- mean(avg_profile[(ceiling(length(avg_profile) / 2) + 1):length(avg_profile)], na.rm = TRUE)
      trend <- if (second_half - first_half > 0.3) {
        "increasing"
      } else if (second_half - first_half < -0.3) {
        "decreasing"
      } else {
        "stable/oscillatory"
      }
      cluster_lines <- c(cluster_lines, paste0("- **Temporal trend:** ", trend))
    }

    ## Format time profile
    profile_str <- paste(
      sprintf("%s: %s", time_points, omicsai::format_num(avg_profile, 2)),
      collapse = ", "
    )
    cluster_lines <- c(cluster_lines, paste0("- **Mean expression profile (scaled):** ", profile_str))

    ## Expression range
    sd_vals <- matrixStats::rowSds(module_expr, na.rm = TRUE)
    cluster_lines <- c(cluster_lines, paste0(
      "- **Mean SD across time:** ", omicsai::format_num(mean(sd_vals, na.rm = TRUE), 2)
    ))
  }

  ## Add top genes by SD
  if (!is.null(data$X) && length(module_genes) > 0) {
    module_X <- data$X[module_genes, , drop = FALSE]
    gene_sd <- matrixStats::rowSds(module_X, na.rm = TRUE)
    names(gene_sd) <- module_genes
    gene_sd <- sort(gene_sd, decreasing = TRUE)
    top_genes <- head(names(gene_sd), 15)

    ## Convert to symbols if possible
    symbols <- top_genes
    if (!is.null(pgx$genes)) {
      symbols <- tryCatch(
        playbase::probe2symbol(top_genes, pgx$genes, "symbol"),
        error = function(e) top_genes
      )
      symbols <- ifelse(is.na(symbols) | symbols == "", top_genes, symbols)
    }
    cluster_lines <- c(
      cluster_lines,
      paste0("- **Top genes (by variability):** ", paste(symbols, collapse = ", "))
    )
  }

  cluster_info <- paste(cluster_lines, collapse = "\n")

  ## -----------------------------------------------------------
  ## Enrichment results (gene set correlations with module)
  ## -----------------------------------------------------------
  ss <- ""
  gset_rho <- data$gset.rho
  if (!is.null(gset_rho) && module %in% colnames(gset_rho)) {
    rho <- gset_rho[, module]
    names(rho) <- rownames(gset_rho)

    ## Compute p-values for correlations
    n_samples <- ncol(pgx$X)
    pv <- tryCatch(
      playbase::cor.pvalue(rho, n_samples),
      error = function(e) rep(NA_real_, length(rho))
    )
    names(pv) <- names(rho)

    ## Select top positively and negatively correlated sets
    ord_pos <- order(rho, decreasing = TRUE)
    ord_neg <- order(rho, decreasing = FALSE)
    n_each <- ceiling(ntop / 2)
    sel <- unique(c(head(ord_pos, n_each), head(ord_neg, n_each)))
    sel <- sel[order(rho[sel], decreasing = TRUE)]

    ## Filter to at least marginally significant
    sig <- sel[!is.na(pv[sel]) & pv[sel] < 0.10]
    if (length(sig) > 0) sel <- sig

    top_gsets <- head(sel, ntop)
    if (length(top_gsets) > 0) {
      pathways <- sub(".*:", "", names(rho)[top_gsets])
      gset_table <- paste0(
        "| Gene Set | Correlation (rho) | P-value |\n",
        "|----------|-------------------|---------|\n",
        paste(
          sprintf(
            "| %s | %s | %s |",
            pathways,
            omicsai::format_num(rho[top_gsets], 3),
            omicsai::format_num(pv[top_gsets], 4)
          ),
          collapse = "\n"
        )
      )

      n_pos <- sum(rho > 0.4 & pv < 0.05, na.rm = TRUE)
      n_neg <- sum(rho < -0.4 & pv < 0.05, na.rm = TRUE)

      ss <- paste0(
        "**Top Correlated Gene Sets:**\n\n",
        gset_table, "\n\n",
        "**Positively correlated (rho > 0.4, p < 0.05):** ", n_pos, " gene sets\n",
        "**Negatively correlated (rho < -0.4, p < 0.05):** ", n_neg, " gene sets"
      )
    }
  }

  if (ss == "") {
    ss <- "No enrichment data available for this module."
  }

  ## -----------------------------------------------------------
  ## Experiment description
  ## -----------------------------------------------------------
  experiment <- ""
  if (!is.null(pgx$description)) {
    experiment <- pgx$description
  }

  list(
    module = module,
    contrast = contrast,
    genesets = ss,
    cluster_info = cluster_info,
    experiment = experiment
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: Time Series AI Summary
# -----------------------------------------------------------------------------

#' Time Series AI Summary UI
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

#' Time Series AI Summary Server
#'
#' Server logic for AI-generated time series cluster summaries.
#' Uses omicsai_text_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param data Reactive returning timeseries_filtered() results
#' @param r_module Reactive returning selected module name
#' @param r_contrast Reactive returning selected contrast name
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_text_server
timeseries_ai_summary_server <- function(id,
                                         pgx,
                                         data,
                                         r_module,
                                         r_contrast,
                                         session,
                                         cache = NULL,
                                         watermark = FALSE) {
  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    res <- data()
    module <- r_module()
    contrast <- r_contrast()
    shiny::req(res, module, contrast)

    timeseries_build_ai_params(
      data = res,
      pgx = pgx,
      module = module,
      contrast = contrast,
      ntop = 20
    )
  })

  # Load template from board.timeseries prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.timeseries/prompts/timeseries_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omicsai::omicsai_load_template(board_template_path)
  })

  # Load TimeSeries methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.timeseries/prompts/TimeSeries_methods.md"
  )
  context_template <- omicsai::omicsai_load_template(context_template_path)

  # Build config reactive from user's selected LLM model
  # Includes board-specific context in system prompt
  config_reactive <- shiny::reactive({
    ai_model <- getUserOption(session, "llm_model")
    shiny::req(ai_model, ai_model != "")

    # Build context params for template substitution
    context_params <- list(
      experiment = pgx$description %||% "omics experiment"
    )

    omicsai::omicsai_config(
      model = ai_model,
      context = context_template,
      context_params = context_params
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
