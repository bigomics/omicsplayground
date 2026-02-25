##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for compare summary
#'
#' Extracts parameters from contrast comparison results to populate the prompt template.
#'
#' @param pgx PGX object
#' @param contrast1 Character; first contrast name
#' @param contrast2 Character; second contrast name
#' @param ntop Integer; number of top genes to include (default 20)
#'
#' @return Named list with template parameters:
#'   contrast1, contrast2, correlation_stats, shared_genes, contrast_specific, experiment
compare_build_ai_params <- function(pgx,
                                    contrast1,
                                    contrast2,
                                    ntop = 20) {

  # Extract experiment description
  experiment <- ""
  if (!is.null(pgx$name) && nzchar(pgx$name)) {
    experiment <- pgx$name
  } else if (!is.null(pgx$description) && nzchar(pgx$description)) {
    experiment <- pgx$description
  }

  # Get fold-change matrices for both contrasts
  meta <- playbase::pgx.getMetaMatrix(pgx)$fc
  fc1 <- NULL
  fc2 <- NULL

  if (contrast1 %in% colnames(meta)) {
    fc1 <- meta[, contrast1, drop = TRUE]
    names(fc1) <- rownames(meta)
  }
  if (contrast2 %in% colnames(meta)) {
    fc2 <- meta[, contrast2, drop = TRUE]
    names(fc2) <- rownames(meta)
  }

  # Map probe IDs to gene symbols if available
  map_to_symbol <- function(probe_names) {
    if (!is.null(pgx$genes) && "symbol" %in% colnames(pgx$genes)) {
      matched <- intersect(probe_names, rownames(pgx$genes))
      if (length(matched) > 0) {
        symbols <- pgx$genes[matched, "symbol"]
        symbols <- ifelse(is.na(symbols) | symbols == "", matched, symbols)
        names(symbols) <- matched
        return(symbols)
      }
    }
    setNames(probe_names, probe_names)
  }

  # Build correlation statistics section
  correlation_stats <- ""
  if (!is.null(fc1) && !is.null(fc2)) {
    gg <- intersect(names(fc1), names(fc2))
    gg <- setdiff(gg, c("", NA, "NA"))

    if (length(gg) > 10) {
      f1 <- fc1[gg]
      f2 <- fc2[gg]

      # Remove NAs for correlation
      valid <- !is.na(f1) & !is.na(f2)
      f1_valid <- f1[valid]
      f2_valid <- f2[valid]

      pearson_r <- cor(f1_valid, f2_valid, method = "pearson")
      spearman_r <- cor(f1_valid, f2_valid, method = "spearman")

      # Count genes significant in each contrast (|logFC| > 0.5 as proxy)
      sig_threshold <- 0.5
      sig1 <- abs(f1_valid) > sig_threshold
      sig2 <- abs(f2_valid) > sig_threshold
      n_sig1 <- sum(sig1)
      n_sig2 <- sum(sig2)
      n_both <- sum(sig1 & sig2)
      n_concordant <- sum(sig1 & sig2 & sign(f1_valid) == sign(f2_valid))
      n_discordant <- sum(sig1 & sig2 & sign(f1_valid) != sign(f2_valid))

      correlation_stats <- paste0(
        "**Fold-change correlation between contrasts:**\n\n",
        "| Metric | Value |\n",
        "|--------|-------|\n",
        "| Pearson r | ", omicsai::omicsai_format_num(pearson_r, 3), " |\n",
        "| Spearman rho | ", omicsai::omicsai_format_num(spearman_r, 3), " |\n",
        "| Total shared genes | ", length(f1_valid), " |\n",
        "| Significant in contrast 1 (|logFC| > ", sig_threshold, ") | ", n_sig1, " |\n",
        "| Significant in contrast 2 (|logFC| > ", sig_threshold, ") | ", n_sig2, " |\n",
        "| Significant in both | ", n_both, " |\n",
        "| Concordant direction | ", n_concordant, " |\n",
        "| Discordant direction | ", n_discordant, " |"
      )
    }
  }

  # Build shared significant genes section
  shared_genes <- ""
  if (!is.null(fc1) && !is.null(fc2)) {
    gg <- intersect(names(fc1), names(fc2))
    gg <- setdiff(gg, c("", NA, "NA"))

    if (length(gg) > 0) {
      f1 <- fc1[gg]
      f2 <- fc2[gg]

      valid <- !is.na(f1) & !is.na(f2)
      f1_valid <- f1[valid]
      f2_valid <- f2[valid]
      gg_valid <- gg[valid]

      # Genes significant in both contrasts
      sig_threshold <- 0.5
      both_sig <- abs(f1_valid) > sig_threshold & abs(f2_valid) > sig_threshold
      shared_idx <- which(both_sig)

      if (length(shared_idx) > 0) {
        # Rank by combined absolute fold change
        combined_fc <- abs(f1_valid[shared_idx]) + abs(f2_valid[shared_idx])
        top_idx <- shared_idx[order(-combined_fc)]
        top_idx <- head(top_idx, ntop)

        top_probes <- gg_valid[top_idx]
        sym_map <- map_to_symbol(top_probes)
        symbols <- sym_map[top_probes]

        direction <- ifelse(
          sign(f1_valid[top_idx]) == sign(f2_valid[top_idx]),
          "concordant", "discordant"
        )

        gene_table <- paste0(
          "| Gene | logFC (contrast 1) | logFC (contrast 2) | Direction |\n",
          "|------|--------------------|--------------------|-----------|\n",
          paste(
            sprintf(
              "| %s | %s | %s | %s |",
              symbols,
              omicsai::omicsai_format_num(f1_valid[top_idx], 3),
              omicsai::omicsai_format_num(f2_valid[top_idx], 3),
              direction
            ),
            collapse = "\n"
          )
        )

        n_concordant <- sum(direction == "concordant")
        n_discordant <- sum(direction == "discordant")

        shared_genes <- paste0(
          "**Top shared significant genes (ranked by combined |logFC|):**\n\n",
          gene_table, "\n\n",
          "**Top ", length(top_idx), " shown:** ",
          n_concordant, " concordant, ", n_discordant, " discordant  \n",
          "**Total shared significant genes:** ", length(shared_idx)
        )
      }
    }
  }

  # Build contrast-specific genes section
  contrast_specific <- ""
  if (!is.null(fc1) && !is.null(fc2)) {
    gg <- intersect(names(fc1), names(fc2))
    gg <- setdiff(gg, c("", NA, "NA"))

    if (length(gg) > 0) {
      f1 <- fc1[gg]
      f2 <- fc2[gg]

      valid <- !is.na(f1) & !is.na(f2)
      f1_valid <- f1[valid]
      f2_valid <- f2[valid]
      gg_valid <- gg[valid]

      sig_threshold <- 0.5
      only1 <- abs(f1_valid) > sig_threshold & abs(f2_valid) <= sig_threshold
      only2 <- abs(f2_valid) > sig_threshold & abs(f1_valid) <= sig_threshold

      ntop_half <- max(floor(ntop / 2), 5)

      parts <- c()

      # Genes unique to contrast 1
      idx1 <- which(only1)
      if (length(idx1) > 0) {
        top1 <- idx1[order(-abs(f1_valid[idx1]))]
        top1 <- head(top1, ntop_half)
        probes1 <- gg_valid[top1]
        sym1 <- map_to_symbol(probes1)[probes1]

        tbl1 <- paste0(
          "**Specific to contrast 1** (", contrast1, "):\n\n",
          "| Gene | logFC (contrast 1) | logFC (contrast 2) |\n",
          "|------|--------------------|--------------------|",
          "\n",
          paste(
            sprintf(
              "| %s | %s | %s |",
              sym1,
              omicsai::omicsai_format_num(f1_valid[top1], 3),
              omicsai::omicsai_format_num(f2_valid[top1], 3)
            ),
            collapse = "\n"
          ),
          "\n\n",
          "**", length(idx1), " genes** significant only in contrast 1 (top ", length(top1), " shown)"
        )
        parts <- c(parts, tbl1)
      }

      # Genes unique to contrast 2
      idx2 <- which(only2)
      if (length(idx2) > 0) {
        top2 <- idx2[order(-abs(f2_valid[idx2]))]
        top2 <- head(top2, ntop_half)
        probes2 <- gg_valid[top2]
        sym2 <- map_to_symbol(probes2)[probes2]

        tbl2 <- paste0(
          "**Specific to contrast 2** (", contrast2, "):\n\n",
          "| Gene | logFC (contrast 1) | logFC (contrast 2) |\n",
          "|------|--------------------|--------------------|",
          "\n",
          paste(
            sprintf(
              "| %s | %s | %s |",
              sym2,
              omicsai::omicsai_format_num(f1_valid[top2], 3),
              omicsai::omicsai_format_num(f2_valid[top2], 3)
            ),
            collapse = "\n"
          ),
          "\n\n",
          "**", length(idx2), " genes** significant only in contrast 2 (top ", length(top2), " shown)"
        )
        parts <- c(parts, tbl2)
      }

      if (length(parts) > 0) {
        contrast_specific <- paste(parts, collapse = "\n\n")
      }
    }
  }

  list(
    contrast1 = contrast1,
    contrast2 = contrast2,
    correlation_stats = correlation_stats,
    shared_genes = shared_genes,
    contrast_specific = contrast_specific,
    experiment = experiment
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: Compare AI Summary
# -----------------------------------------------------------------------------


#' Compare AI Summary Server
#'
#' Server logic for compare AI-generated summaries.
#' Uses omicsai_text_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param contrast1 Reactive returning selected first contrast
#' @param contrast2 Reactive returning selected second contrast
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_text_server
compare_ai_summary_server <- function(id,
                                      pgx,
                                      contrast1,
                                      contrast2,
                                      session,
                                      cache = NULL,
                                      watermark = FALSE) {

  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    ct1 <- contrast1()
    ct2 <- contrast2()
    shiny::req(ct1, ct2)

    compare_build_ai_params(
      pgx = pgx,
      contrast1 = ct1,
      contrast2 = ct2,
      ntop = 20
    )
  })

  # Load template from board.compare prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.compare/prompts/compare_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omicsai::omicsai_load_template(board_template_path)
  })

  # Load compare methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.compare/prompts/Compare_methods.md"
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
