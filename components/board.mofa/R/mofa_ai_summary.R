##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for MOFA factor summary
#'
#' Extracts parameters from MOFA results to populate the prompt template.
#' Collects enrichment data, top loading features, and trait associations
#' for the selected factor.
#'
#' @param mofa MOFA results object
#' @param pgx PGX object (for gene annotation and experiment info)
#' @param factor_name Character; factor name (e.g., "Factor1")
#' @param ntop Integer; number of top genes/sets to include (default 20)
#'
#' @return Named list with template parameters:
#'   factor, genesets, top_genes, traits, experiment
mofa_build_ai_params <- function(mofa,
                                 pgx,
                                 factor_name,
                                 ntop = 20) {

  # --- Extract enriched pathways ---
  genesets <- ""
  if (!is.null(mofa$gsea) && !is.null(mofa$gsea$table) &&
    factor_name %in% names(mofa$gsea$table)) {
    df <- mofa$gsea$table[[factor_name]]
    if (is.data.frame(df) && nrow(df) > 0 &&
      all(c("pathway", "NES", "padj") %in% colnames(df))) {
      df <- df[order(-abs(df$NES)), , drop = FALSE]
      sig <- df[!is.na(df$padj) & df$padj < 0.05, , drop = FALSE]
      top_gse <- head(sig, 8)
      if (nrow(top_gse) > 0) {
        pathways <- sub(".*:", "", top_gse$pathway)

        # Build size column (may have multiple size columns for multi-omics)
        size_cols <- grep("^size", colnames(top_gse), value = TRUE)
        if (length(size_cols) > 0) {
          sizes <- apply(top_gse[, size_cols, drop = FALSE], 1, function(x) {
            paste(x, collapse = "/")
          })
        } else {
          sizes <- rep("NA", nrow(top_gse))
        }

        gset_table <- paste0(
          "| Pathway | NES | padj | Size |\n",
          "|---------|-----|------|------|\n",
          paste(
            sprintf(
              "| %s | %s | %s | %s |",
              pathways,
              omicsai::format_num(top_gse$NES, 2),
              omicsai::format_num(top_gse$padj, 3),
              sizes
            ),
            collapse = "\n"
          )
        )

        n_sig <- sum(df$padj < 0.05, na.rm = TRUE)
        genesets <- paste0(
          "**Top Enriched Pathways (padj < 0.05):**\n\n",
          gset_table, "\n\n",
          "**Total significant pathways:** ", n_sig, " of ", nrow(df), " tested"
        )
      }
    }
  }

  # --- Extract top loading genes/features ---
  top_genes <- ""
  if (!is.null(mofa$W) && factor_name %in% colnames(mofa$W)) {
    w <- mofa$W[, factor_name]
    w <- w[!is.na(w)]
    w <- sort(abs(w), decreasing = TRUE)
    top_features <- head(names(w), ntop)
    top_weights <- head(w, ntop)

    # Try to get gene symbols from annotation
    symbols <- top_features
    if (!is.null(pgx$genes)) {
      symbols <- tryCatch(
        playbase::probe2symbol(top_features, pgx$genes, "symbol"),
        error = function(e) top_features
      )
      symbols <- ifelse(is.na(symbols) | symbols == "", top_features, symbols)
    }

    # Try to get centrality data
    centrality_vals <- rep(NA_real_, length(top_features))
    centrality_data <- tryCatch(
      playbase::mofa.plot_centrality(mofa, factor_name, justdata = TRUE),
      error = function(e) NULL
    )
    if (is.data.frame(centrality_data) && "centrality" %in% colnames(centrality_data) &&
      "feature" %in% colnames(centrality_data)) {
      m <- match(top_features, centrality_data$feature)
      centrality_vals <- centrality_data$centrality[m]
    }

    gene_table <- paste0(
      "| Feature | Symbol | Weight | Centrality |\n",
      "|---------|--------|--------|------------|\n",
      paste(
        sprintf(
          "| %s | %s | %s | %s |",
          sub(";.*", ";[...]", top_features),
          symbols,
          omicsai::format_num(top_weights, 3),
          omicsai::format_num(centrality_vals, 2)
        ),
        collapse = "\n"
      )
    )

    top_genes <- paste0(
      "**Top Loading Features (ranked by absolute weight):**\n\n",
      gene_table, "\n\n",
      "**Top ", length(top_features), " of ",
      length(w), " features shown**"
    )
  }

  # --- Extract associated traits ---
  traits <- "None"
  if (!is.null(mofa$Y) && !is.null(mofa$F) &&
    factor_name %in% colnames(mofa$F)) {
    # Compute factor-trait correlations
    factor_scores <- mofa$F[, factor_name]
    trait_names <- colnames(mofa$Y)
    if (length(trait_names) > 0) {
      cors <- sapply(trait_names, function(t) {
        y <- mofa$Y[, t]
        if (is.numeric(y) && sum(!is.na(y)) > 3) {
          tryCatch(
            cor(factor_scores, y, use = "pairwise.complete.obs"),
            error = function(e) NA_real_
          )
        } else {
          NA_real_
        }
      })
      cors <- cors[!is.na(cors)]
      if (length(cors) > 0) {
        # Report traits with |cor| > 0.3, sorted by absolute correlation
        sig_traits <- cors[abs(cors) > 0.3]
        if (length(sig_traits) > 0) {
          sig_traits <- sort(abs(sig_traits), decreasing = TRUE)
          traits <- paste(
            sprintf("%s (r=%s)", names(sig_traits), omicsai::format_num(sig_traits, 2)),
            collapse = ", "
          )
        } else {
          # If no strong correlations, show top 3
          top_cors <- head(sort(abs(cors), decreasing = TRUE), 3)
          traits <- paste(
            sprintf("%s (r=%s)", names(top_cors), omicsai::format_num(top_cors, 2)),
            collapse = ", "
          )
        }
      }
    }
  }
  # Fallback: use sample phenotype column names
  if (traits == "None" && !is.null(mofa$samples)) {
    traits <- paste(head(colnames(mofa$samples), 5), collapse = ", ")
  }

  # --- Get experiment description ---
  experiment <- mofa$experiment %||% ""
  if (experiment == "" && !is.null(pgx$annotation)) {
    experiment <- pgx$annotation %||% ""
  }

  list(
    factor = factor_name,
    genesets = genesets,
    top_genes = top_genes,
    traits = traits,
    experiment = experiment
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: MOFA AI Summary
# -----------------------------------------------------------------------------

#' MOFA AI Summary UI
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

#' MOFA AI Summary Server
#'
#' Server logic for MOFA AI-generated factor summaries.
#' Uses omicsai_text_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param mofa Reactive returning MOFA results object
#' @param pgx PGX object (reactive values)
#' @param r_factor Reactive returning selected factor name
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_text_server
mofa_ai_summary_server <- function(id,
                                   mofa,
                                   pgx,
                                   r_factor,
                                   session,
                                   cache = NULL,
                                   watermark = FALSE) {
  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    m <- mofa()
    factor_name <- r_factor()
    shiny::req(m, factor_name)

    mofa_build_ai_params(
      mofa = m,
      pgx = pgx,
      factor_name = factor_name,
      ntop = 20
    )
  })

  # Load template from board.mofa prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.mofa/prompts/mofa_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omicsai::omicsai_load_template(board_template_path)
  })

  # Load MOFA methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.mofa/prompts/MOFA_methods.md"
  )
  context_template <- omicsai::omicsai_load_template(context_template_path)

  # Build config reactive from user's selected LLM model
  # Includes board-specific context in system prompt
  config_reactive <- shiny::reactive({
    ai_model <- getUserOption(session, "llm_model")
    m <- mofa()
    shiny::req(ai_model, ai_model != "", m)

    # Build context params for template substitution
    context_params <- list(
      experiment = m$experiment %||% "multi-omics experiment"
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
