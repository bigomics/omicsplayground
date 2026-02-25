##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for feature map summary
#'
#' Extracts parameters from feature map data to populate the prompt template.
#'
#' @param pgx PGX object
#' @param contrast Character; contrast name (e.g., "groupA_vs_groupB")
#' @param feature_level Character; "gene" or "geneset"
#' @param ntop Integer; number of top features to include (default 20)
#'
#' @return Named list with template parameters:
#'   contrast, feature_level, top_features, nearby_features, experiment
featuremap_build_ai_params <- function(pgx,
                                       contrast,
                                       feature_level = "gene",
                                       ntop = 20) {

  # Extract experiment description
  experiment <- ""
  if (!is.null(pgx$name) && nzchar(pgx$name)) {
    experiment <- pgx$name
  } else if (!is.null(pgx$description) && nzchar(pgx$description)) {
    experiment <- pgx$description
  }

  # Select data based on feature level
  if (feature_level == "geneset") {
    level <- "geneset"
    pos <- pgx$cluster.gsets$pos[["umap2d"]]
  } else {
    level <- "gene"
    pos <- pgx$cluster.genes$pos[["umap2d"]]
  }

  # Get fold change matrix
  F_mat <- playbase::pgx.getMetaMatrix(pgx, level = level)$fc
  F_scaled <- scale(F_mat, center = FALSE)
  rms_fc <- sqrt(rowMeans(F_scaled**2, na.rm = TRUE))

  # Conform to UMAP positions
  gg <- intersect(rownames(pos), names(rms_fc))
  pos <- pos[gg, , drop = FALSE]
  rms_fc <- rms_fc[gg]
  F_mat <- F_mat[gg, , drop = FALSE]

  # Get contrast-specific fold changes
  contrast_fc <- NULL
  if (contrast %in% colnames(F_mat)) {
    contrast_fc <- F_mat[, contrast]
  }

  # Build top features table (ranked by rms.FC)
  top_features <- ""
  ord <- order(-rms_fc)
  top_idx <- head(ord, ntop)
  top_names <- names(rms_fc)[top_idx]

  if (length(top_names) > 0) {
    # Get gene symbols for gene-level features
    display_names <- top_names
    if (feature_level == "gene" && !is.null(pgx$genes) && "symbol" %in% colnames(pgx$genes)) {
      matched <- intersect(top_names, rownames(pgx$genes))
      if (length(matched) > 0) {
        symbols <- pgx$genes[matched, "symbol"]
        symbols <- ifelse(is.na(symbols) | symbols == "", matched, symbols)
        display_names[top_names %in% matched] <- symbols
      }
    } else if (feature_level == "geneset") {
      # Strip collection prefix for readability
      display_names <- sub(".*:", "", top_names)
    }

    # Build contrast FC column if available
    if (!is.null(contrast_fc)) {
      fc_vals <- contrast_fc[top_names]
      feat_table <- paste0(
        "| Feature | rms.FC | logFC (", contrast, ") |\n",
        "|---------|--------|", paste0(rep("-", nchar(contrast) + 9), collapse = ""), "|\n",
        paste(
          sprintf(
            "| %s | %s | %s |",
            display_names,
            omicsai::omicsai_format_num(rms_fc[top_names], 3),
            omicsai::omicsai_format_num(fc_vals, 3)
          ),
          collapse = "\n"
        )
      )
    } else {
      feat_table <- paste0(
        "| Feature | rms.FC |\n",
        "|---------|--------|\n",
        paste(
          sprintf(
            "| %s | %s |",
            display_names,
            omicsai::omicsai_format_num(rms_fc[top_names], 3)
          ),
          collapse = "\n"
        )
      )
    }

    n_total <- length(rms_fc)
    top_features <- paste0(
      "**Top Variable Features (by rms.FC):**\n\n",
      feat_table, "\n\n",
      "**Total features mapped:** ", n_total
    )
  }

  # Build nearby features section
  # For the top 5 features, find their nearest neighbors on the UMAP
  nearby_features <- ""
  if (nrow(pos) > 0 && length(top_names) > 0) {
    anchor_names <- head(top_names, 5)
    nearby_lines <- c()

    for (anchor in anchor_names) {
      if (!anchor %in% rownames(pos)) next
      anchor_pos <- pos[anchor, , drop = FALSE]
      dists <- sqrt(rowSums((pos - matrix(anchor_pos, nrow = nrow(pos), ncol = 2, byrow = TRUE))^2))
      names(dists) <- rownames(pos)
      dists <- dists[names(dists) != anchor]
      dists <- sort(dists)
      neighbors <- head(names(dists), 10)

      # Get display names for neighbors
      if (feature_level == "gene" && !is.null(pgx$genes) && "symbol" %in% colnames(pgx$genes)) {
        matched <- intersect(neighbors, rownames(pgx$genes))
        neigh_display <- neighbors
        if (length(matched) > 0) {
          syms <- pgx$genes[matched, "symbol"]
          syms <- ifelse(is.na(syms) | syms == "", matched, syms)
          neigh_display[neighbors %in% matched] <- syms
        }
        anchor_display <- anchor
        if (anchor %in% rownames(pgx$genes)) {
          s <- pgx$genes[anchor, "symbol"]
          if (!is.na(s) && nzchar(s)) anchor_display <- s
        }
      } else if (feature_level == "geneset") {
        neigh_display <- sub(".*:", "", neighbors)
        anchor_display <- sub(".*:", "", anchor)
      } else {
        neigh_display <- neighbors
        anchor_display <- anchor
      }

      # Get fold changes for neighbors
      neigh_fc <- ""
      if (!is.null(contrast_fc)) {
        fc_vals <- contrast_fc[neighbors]
        neigh_fc <- paste0(" (logFC: ", paste(omicsai::omicsai_format_num(fc_vals, 3), collapse = ", "), ")")
      }

      nearby_lines <- c(
        nearby_lines,
        paste0(
          "- **", anchor_display, "** neighbors: ",
          paste(neigh_display, collapse = ", "),
          neigh_fc
        )
      )
    }

    if (length(nearby_lines) > 0) {
      nearby_features <- paste0(
        "**Nearest neighbors on the UMAP for the top variable features:**\n\n",
        paste(nearby_lines, collapse = "\n")
      )
    }
  }

  list(
    contrast = contrast,
    feature_level = feature_level,
    top_features = top_features,
    nearby_features = nearby_features,
    experiment = experiment
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: Feature Map AI Summary
# -----------------------------------------------------------------------------

#' Feature Map AI Summary UI
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

#' Feature Map AI Summary Server
#'
#' Server logic for feature map AI-generated summaries.
#' Uses omicsai_text_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param sigvar Reactive returning selected comparison(s) or phenotype variable(s)
#' @param feature_level Reactive returning "gene" or "geneset"
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_text_server
featuremap_ai_summary_server <- function(id,
                                         pgx,
                                         sigvar,
                                         feature_level,
                                         session,
                                         cache = NULL,
                                         watermark = FALSE) {

  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    sv <- sigvar()
    fl <- feature_level()
    shiny::req(sv, fl)

    # Use first comparison/phenotype as the contrast label
    contrast <- sv[1]

    featuremap_build_ai_params(
      pgx = pgx,
      contrast = contrast,
      feature_level = fl,
      ntop = 20
    )
  })

  # Load template from board.featuremap prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.featuremap/prompts/featuremap_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omicsai::omicsai_load_template(board_template_path)
  })

  # Load feature map methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.featuremap/prompts/FeatureMap_methods.md"
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
