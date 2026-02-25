##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for TCGA survival summary
#'
#' Extracts parameters from the PGX object and contrast to populate
#' the TCGA summary prompt template.
#'
#' @param pgx PGX object
#' @param contrast Character; contrast name (e.g., "groupA_vs_groupB")
#' @param sigtype Character; signature type ("contrast" or "genelist")
#' @param genelist Character; pasted gene list (used when sigtype == "genelist")
#' @param ntop Integer; number of top genes to include (default 15)
#'
#' @return Named list with template parameters:
#'   contrast, phenotype, signature_type, top_genes, survival_results, experiment
tcga_build_ai_params <- function(pgx,
                                 contrast,
                                 sigtype = "contrast",
                                 genelist = NULL,
                                 ntop = 15) {

  # Extract experiment description
  experiment <- pgx$name %||% pgx$description %||% "omics experiment"

  # Extract phenotype from contrast name (readable form)
  phenotype <- gsub("_", " ", contrast)

  # Determine signature type label
  signature_type <- if (sigtype == "contrast") {
    paste0("Contrast-based (", contrast, ")")
  } else {
    "Custom gene list"
  }

  # Build top genes section
  top_genes <- ""
  if (sigtype == "contrast") {
    # Get fold-change signature from contrast
    res <- tryCatch(
      playbase::pgx.getMetaFoldChangeMatrix(pgx, what = "meta"),
      error = function(e) NULL
    )

    if (!is.null(res) && contrast %in% colnames(res$fc)) {
      sig <- res$fc[, contrast]
      sig <- sig[!is.na(sig)]
      sig <- sort(sig, decreasing = TRUE)

      # Get top upregulated and downregulated genes
      top_up <- head(sig[sig > 0], ntop)
      top_down <- tail(sig[sig < 0], ntop)
      top_all <- c(top_up, rev(top_down))

      if (length(top_all) > 0) {
        # Map probe names to gene symbols if possible
        probe_names <- names(top_all)
        symbols <- probe_names
        if (!is.null(pgx$genes) && "symbol" %in% colnames(pgx$genes)) {
          matched <- intersect(probe_names, rownames(pgx$genes))
          if (length(matched) > 0) {
            sym <- pgx$genes[matched, "symbol"]
            symbols[match(matched, probe_names)] <- ifelse(
              is.na(sym) | sym == "", matched, sym
            )
          }
        }

        # Deduplicate by symbol, keeping highest absolute FC
        df <- data.frame(
          symbol = symbols,
          fc = as.numeric(top_all),
          stringsAsFactors = FALSE
        )
        df <- df[order(-abs(df$fc)), , drop = FALSE]
        df <- df[!duplicated(df$symbol), , drop = FALSE]
        df <- head(df, ntop)

        if (nrow(df) > 0) {
          n_up <- sum(df$fc > 0, na.rm = TRUE)
          n_down <- sum(df$fc < 0, na.rm = TRUE)

          gene_table <- paste0(
            "| Gene | logFC |\n",
            "|------|-------|\n",
            paste(
              sprintf("| %s | %s |", df$symbol, omicsai::omicsai_format_num(df$fc, 3)),
              collapse = "\n"
            )
          )

          top_genes <- paste0(
            "**Top signature genes (ranked by |logFC|):**\n\n",
            gene_table, "\n\n",
            "**Direction:** ", n_up, " upregulated, ", n_down, " downregulated (top ", nrow(df), " shown)"
          )
        }
      }
    }
  } else if (sigtype == "genelist" && !is.null(genelist) && nzchar(genelist)) {
    genes <- as.character(genelist)
    genes <- strsplit(genes, split = "[\t, \n]")[[1]]
    genes <- gsub("[ ]", "", genes)
    genes <- genes[nzchar(genes)]

    if (length(genes) > 0) {
      genes_shown <- head(genes, 30)
      top_genes <- paste0(
        "**Custom gene list (", length(genes), " genes):**\n\n",
        paste(genes_shown, collapse = ", "),
        if (length(genes) > 30) paste0("\n... and ", length(genes) - 30, " more") else ""
      )
    }
  }

  # Build survival results section
  # Note: The actual survival p-values come from pgx.testTCGAsurvival which
  # generates a plot. We describe the analysis setup since the computed
  # survival results are rendered as a plot (not returned as data).
  survival_results <- ""

  # Try to extract available TCGA-related metadata if stored in pgx
  n_genes_in_sig <- 0
  if (sigtype == "contrast") {
    res <- tryCatch(
      playbase::pgx.getMetaFoldChangeMatrix(pgx, what = "meta"),
      error = function(e) NULL
    )
    if (!is.null(res) && contrast %in% colnames(res$fc)) {
      sig <- res$fc[, contrast]
      n_genes_in_sig <- sum(!is.na(sig))
    }
  } else if (!is.null(genelist) && nzchar(genelist)) {
    genes <- strsplit(as.character(genelist), split = "[\t, \n]")[[1]]
    genes <- gsub("[ ]", "", genes)
    n_genes_in_sig <- sum(nzchar(genes))
  }

  survival_results <- paste0(
    "The TCGA survival analysis was performed using the Kaplan-Meier method across 32 TCGA cancer types.\n\n",
    "**Signature size:** ", n_genes_in_sig, " genes in the input signature\n\n",
    "Each TCGA cohort (e.g., BRCA, LUAD, COAD, GBM, etc.) was dichotomized into patients ",
    "positively and negatively correlated with the input gene signature. ",
    "Survival probabilities were computed and tested using the log-rank test.\n\n",
    "Please interpret the biological relevance of this gene signature in the context of ",
    "cancer survival, considering which cancer types are most likely to show significant ",
    "survival associations based on the signature genes and their functions."
  )

  list(
    contrast = contrast,
    phenotype = phenotype,
    signature_type = signature_type,
    top_genes = top_genes,
    survival_results = survival_results,
    experiment = experiment
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: TCGA AI Summary
# -----------------------------------------------------------------------------

#' TCGA AI Summary UI
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

#' TCGA AI Summary Server
#'
#' Server logic for TCGA AI-generated summaries.
#' Uses omicsai_text_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param contrast Reactive returning selected contrast
#' @param sigtype Reactive returning signature type ("contrast" or "genelist")
#' @param genelist Reactive returning pasted gene list
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_text_server
tcga_ai_summary_server <- function(id,
                                   pgx,
                                   contrast,
                                   sigtype,
                                   genelist,
                                   session,
                                   cache = NULL,
                                   watermark = FALSE) {

  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    ct <- contrast()
    st <- sigtype()
    gl <- genelist()
    shiny::req(ct)

    tcga_build_ai_params(
      pgx = pgx,
      contrast = ct,
      sigtype = st,
      genelist = gl,
      ntop = 15
    )
  })

  # Load template from board.tcga prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.tcga/prompts/tcga_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omicsai::omicsai_load_template(board_template_path)
  })

  # Load TCGA methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.tcga/prompts/TCGA_methods.md"
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
