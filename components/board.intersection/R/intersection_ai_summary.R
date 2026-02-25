##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for intersection summary
#'
#' Extracts parameters from intersection analysis results to populate the
#' prompt template.
#'
#' @param pgx PGX object
#' @param comparisons Character vector; selected contrast names
#' @param level Character; feature level ("gene" or "geneset")
#' @param fc_matrix List with elements fc, qv (from getFoldChangeMatrix)
#' @param active_fc_matrix List with elements fc, qv (from getActiveFoldChangeMatrix)
#' @param fdr Numeric; FDR threshold used for significance calls
#' @param lfc Numeric; logFC threshold used for significance calls
#' @param ntop Integer; number of top genes to include in heatmap summary (default 30)
#'
#' @return Named list with template parameters:
#'   contrasts, level, venn_summary, heatmap_genes, correlation_info, n_contrasts, experiment
intersection_build_ai_params <- function(pgx,
                                         comparisons,
                                         level,
                                         fc_matrix,
                                         active_fc_matrix,
                                         fdr = 0.05,
                                         lfc = 0.2,
                                         ntop = 30) {

  # Extract experiment description
  experiment <- pgx$name %||% pgx$description %||% "omics experiment"

  # Format selected contrasts
  contrasts_str <- paste(comparisons, collapse = ", ")
  n_contrasts <- length(comparisons)

  # --- Venn diagram summary ---
  venn_summary <- ""
  if (!is.null(fc_matrix) && !is.null(fc_matrix$fc) && !is.null(fc_matrix$qv)) {
    fc <- fc_matrix$fc
    qv <- fc_matrix$qv

    # Subset to selected comparisons
    sel <- intersect(comparisons, colnames(fc))
    if (length(sel) >= 2) {
      fc_sel <- fc[, sel, drop = FALSE]
      qv_sel <- qv[, sel, drop = FALSE]

      # Compute significance calls
      sig <- (qv_sel <= fdr & abs(fc_sel) >= lfc)
      sig[is.na(sig)] <- FALSE

      # Compute overlap counts
      n_sig_per_contrast <- colSums(sig)
      n_any <- sum(rowSums(sig) >= 1)
      n_all <- sum(rowSums(sig) == ncol(sig))

      # Per-contrast significance counts
      contrast_counts <- paste(
        sprintf("- **%s**: %d significant features", sel, n_sig_per_contrast),
        collapse = "\n"
      )

      # Pairwise overlaps
      pairwise_lines <- c()
      if (length(sel) >= 2) {
        pairs <- utils::combn(sel, 2, simplify = FALSE)
        for (p in pairs) {
          n_shared <- sum(sig[, p[1]] & sig[, p[2]])
          pairwise_lines <- c(pairwise_lines, sprintf(
            "- **%s** & **%s**: %d shared significant features",
            p[1], p[2], n_shared
          ))
        }
      }
      pairwise_str <- paste(pairwise_lines, collapse = "\n")

      venn_summary <- paste0(
        "**Significance thresholds:** FDR < ", fdr, ", |logFC| >= ", lfc, "\n\n",
        "**Per-contrast counts:**\n", contrast_counts, "\n\n",
        "**Pairwise overlaps:**\n", pairwise_str, "\n\n",
        "**Features significant in at least one contrast:** ", n_any, "\n",
        "**Features significant in ALL ", n_contrasts, " contrasts:** ", n_all
      )
    }
  }

  # --- Heatmap top genes ---
  heatmap_genes <- ""
  if (!is.null(active_fc_matrix) && !is.null(active_fc_matrix$fc)) {
    FC <- active_fc_matrix$fc
    if (ncol(FC) > 0 && nrow(FC) > 0) {
      # Order by mean squared FC (most variable)
      FC <- FC[order(-rowMeans(FC^2, na.rm = TRUE)), , drop = FALSE]
      FC <- head(FC, ntop)

      # Convert probe IDs to gene symbols if at gene level
      row_labels <- rownames(FC)
      if (level == "gene" && !is.null(pgx$genes)) {
        matched <- intersect(rownames(FC), rownames(pgx$genes))
        if (length(matched) > 0) {
          symbols <- pgx$genes[matched, "gene_name"]
          valid <- !is.na(symbols) & nzchar(symbols)
          row_labels[match(matched, rownames(FC))] <- ifelse(
            valid, symbols, matched
          )
        }
      }

      # Build markdown table
      header <- paste0(
        "| Feature | ", paste(colnames(FC), collapse = " | "), " |\n",
        "|---------|", paste(rep("------", ncol(FC)), collapse = "|"), "|\n"
      )
      rows <- vapply(seq_len(nrow(FC)), function(i) {
        vals <- paste(omicsai::omicsai_format_num(FC[i, ], 3), collapse = " | ")
        sprintf("| %s | %s |", row_labels[i], vals)
      }, character(1))

      heatmap_genes <- paste0(
        "**Top ", nrow(FC), " features ranked by mean squared fold change:**\n\n",
        header,
        paste(rows, collapse = "\n")
      )
    }
  }

  # --- Correlation info ---
  correlation_info <- ""
  if (!is.null(fc_matrix) && !is.null(fc_matrix$fc)) {
    fc_all <- fc_matrix$fc
    sel <- intersect(comparisons, colnames(fc_all))
    if (length(sel) >= 2) {
      fc_sub <- fc_all[, sel, drop = FALSE]
      fc_sub[is.na(fc_sub)] <- 0
      # Use top 1000 most variable features
      jj <- head(order(-rowMeans(fc_sub^2, na.rm = TRUE)), 1000)
      fc_ranked <- apply(fc_sub[jj, , drop = FALSE], 2, rank, na.last = "keep")
      R <- cor(fc_ranked, use = "pairwise")
      R <- round(R, digits = 2)

      # Build correlation table
      header <- paste0(
        "| Contrast | ", paste(colnames(R), collapse = " | "), " |\n",
        "|----------|", paste(rep("------", ncol(R)), collapse = "|"), "|\n"
      )
      rows <- vapply(seq_len(nrow(R)), function(i) {
        vals <- paste(omicsai::omicsai_format_num(R[i, ], 2), collapse = " | ")
        sprintf("| %s | %s |", rownames(R)[i], vals)
      }, character(1))

      correlation_info <- paste0(
        "**Pairwise Pearson correlation matrix (rank-based, top 1000 features):**\n\n",
        header,
        paste(rows, collapse = "\n")
      )
    }
  }

  list(
    contrasts = contrasts_str,
    level = level,
    venn_summary = venn_summary,
    heatmap_genes = heatmap_genes,
    correlation_info = correlation_info,
    n_contrasts = as.character(n_contrasts),
    experiment = experiment
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: Intersection AI Summary
# -----------------------------------------------------------------------------

#' Intersection AI Summary UI
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

#' Intersection AI Summary Server
#'
#' Server logic for intersection AI-generated summaries.
#' Uses omicsai_text_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param input_comparisons Reactive returning selected contrast names
#' @param level Reactive returning feature level ("gene" or "geneset")
#' @param getFoldChangeMatrix Reactive returning full fold-change matrix list
#' @param getActiveFoldChangeMatrix Reactive returning active fold-change matrix list
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_text_server
intersection_ai_summary_server <- function(id,
                                           pgx,
                                           input_comparisons,
                                           level,
                                           getFoldChangeMatrix,
                                           getActiveFoldChangeMatrix,
                                           session,
                                           cache = NULL,
                                           watermark = FALSE) {

  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    comparisons <- input_comparisons()
    lvl <- level()
    fc_matrix <- getFoldChangeMatrix()
    active_fc_matrix <- getActiveFoldChangeMatrix()
    shiny::req(comparisons, fc_matrix)

    intersection_build_ai_params(
      pgx = pgx,
      comparisons = comparisons,
      level = lvl,
      fc_matrix = fc_matrix,
      active_fc_matrix = active_fc_matrix,
      fdr = 0.05,
      lfc = 0.2,
      ntop = 30
    )
  })

  # Load template from board.intersection prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.intersection/prompts/intersection_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omicsai::omicsai_load_template(board_template_path)
  })

  # Load intersection methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.intersection/prompts/Intersection_methods.md"
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

  AiTextCardServer(
    id = id,
    params_reactive = params_reactive,
    template_reactive = template_reactive,
    config_reactive = config_reactive,
    cache = cache,
    watermark = watermark
  )
}
