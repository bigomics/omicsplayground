##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for single-cell summary
#'
#' Extracts parameters from single-cell / cell profiling results to populate the
#' prompt template.
#'
#' @param pgx PGX object
#' @param refset Character; name of the reference dataset (e.g., "LM22")
#' @param dcmethod Character; deconvolution method name (e.g., "meta.prod")
#' @param ntop Integer; number of top cell types / marker genes to include (default 20)
#'
#' @return Named list with template parameters:
#'   refset, dcmethod, cell_types, marker_genes, phenotype_info, experiment
singlecell_build_ai_params <- function(pgx,
                                       refset,
                                       dcmethod,
                                       ntop = 20) {

  # Extract experiment description
  experiment <- pgx$name %||% pgx$description %||% "omics experiment"

  # ---- Cell type proportions ----
  cell_types <- ""
  deconv_scores <- NULL
  if ("deconv" %in% names(pgx) &&
    refset %in% names(pgx$deconv) &&
    dcmethod %in% names(pgx$deconv[[refset]])) {
    deconv_scores <- pgx$deconv[[refset]][[dcmethod]]
    deconv_scores <- pmax(deconv_scores, 0)

    # Normalize to proportions per sample
    row_sums <- rowSums(deconv_scores, na.rm = TRUE)
    props <- deconv_scores / (1e-20 + row_sums)

    # Mean proportion across all samples
    mean_props <- colMeans(props, na.rm = TRUE)
    mean_props <- sort(mean_props, decreasing = TRUE)
    mean_props <- head(mean_props, ntop)

    # Build per-group breakdown if group info available
    group_breakdown <- ""
    if (!is.null(pgx$model.parameters$group)) {
      grp <- pgx$model.parameters$group[rownames(props)]
      grp_levels <- unique(grp)
      if (length(grp_levels) >= 2 && length(grp_levels) <= 10) {
        top_types <- names(mean_props)
        grp_means <- do.call(rbind, lapply(grp_levels, function(g) {
          idx <- which(grp == g)
          if (length(idx) > 0) colMeans(props[idx, top_types, drop = FALSE], na.rm = TRUE)
          else rep(0, length(top_types))
        }))
        rownames(grp_means) <- grp_levels
        colnames(grp_means) <- top_types

        header <- paste0(
          "| Cell Type | Overall |",
          paste(grp_levels, collapse = " | "), " |\n",
          "|-----------|---------|",
          paste(rep("---", length(grp_levels)), collapse = "|"), "|\n"
        )
        rows <- vapply(seq_along(top_types), function(i) {
          ct <- top_types[i]
          overall <- omicsai::format_num(mean_props[ct] * 100, 1)
          per_grp <- vapply(grp_levels, function(g) {
            omicsai::format_num(grp_means[g, ct] * 100, 1)
          }, character(1))
          paste0("| ", ct, " | ", overall, "% |",
            paste(paste0(per_grp, "%"), collapse = " | "), " |")
        }, character(1))

        group_breakdown <- paste0(
          "\n\n**Cell type proportions by group (%):**\n\n",
          header,
          paste(rows, collapse = "\n")
        )
      }
    }

    # Summary statistics
    n_detected <- sum(mean_props > 0.01)

    cell_types <- paste0(
      "**Top Cell Types (by mean proportion across samples):**\n\n",
      "| Cell Type | Mean Proportion |\n",
      "|-----------|----------------|\n",
      paste(
        sprintf("| %s | %s%% |", names(mean_props), omicsai::format_num(mean_props * 100, 1)),
        collapse = "\n"
      ), "\n\n",
      "**Cell types with >1% mean proportion:** ", n_detected,
      " of ", ncol(deconv_scores), " in reference",
      group_breakdown
    )
  }

  # ---- Marker genes ----
  marker_genes <- ""
  if (!is.null(pgx$X) && nrow(pgx$X) > 0) {
    # Compute top variable genes (same logic as markersplot)
    X <- pgx$X
    gene_sd <- apply(X, 1, sd, na.rm = TRUE)
    top_idx <- head(order(-gene_sd), min(ntop, length(gene_sd)))
    top_probes <- rownames(X)[top_idx]

    # Map to symbols if possible
    symbols <- top_probes
    if (!is.null(pgx$genes) && "symbol" %in% colnames(pgx$genes)) {
      matched <- intersect(top_probes, rownames(pgx$genes))
      if (length(matched) > 0) {
        sym <- pgx$genes[matched, "symbol"]
        sym[is.na(sym) | sym == ""] <- matched[is.na(sym) | sym == ""]
        symbols[match(matched, top_probes)] <- sym
      }
    }

    sd_vals <- gene_sd[top_probes]
    mean_vals <- rowMeans(X[top_probes, , drop = FALSE], na.rm = TRUE)

    marker_table <- paste0(
      "| Gene | Mean Expression | SD |\n",
      "|------|----------------|----|\n",
      paste(
        sprintf("| %s | %s | %s |", symbols, omicsai::format_num(mean_vals, 2), omicsai::format_num(sd_vals, 2)),
        collapse = "\n"
      )
    )

    marker_genes <- paste0(
      "**Top ", length(top_probes), " most variable genes (ranked by SD):**\n\n",
      marker_table, "\n\n",
      "**Total genes in dataset:** ", nrow(X)
    )
  }

  # ---- Phenotype information ----
  phenotype_info <- ""
  if (!is.null(pgx$samples) && ncol(pgx$samples) > 0) {
    pheno_cols <- colnames(pgx$samples)
    pheno_cols <- grep("sample|donor|id|batch", pheno_cols,
      invert = TRUE, value = TRUE, ignore.case = TRUE
    )
    pheno_cols <- head(pheno_cols, 5)

    if (length(pheno_cols) > 0) {
      pheno_summaries <- vapply(pheno_cols, function(col) {
        vals <- pgx$samples[[col]]
        vals <- vals[!is.na(vals) & vals != "" & vals != "NA"]
        n_unique <- length(unique(vals))
        if (n_unique <= 10) {
          tbl <- sort(table(vals), decreasing = TRUE)
          counts <- paste(sprintf("%s (n=%d)", names(tbl), tbl), collapse = ", ")
          sprintf("- **%s**: %s", col, counts)
        } else {
          sprintf("- **%s**: %d unique values", col, n_unique)
        }
      }, character(1))

      n_samples <- nrow(pgx$samples)
      phenotype_info <- paste0(
        "**Sample phenotypes (n=", n_samples, "):**\n\n",
        paste(pheno_summaries, collapse = "\n")
      )
    }
  }

  list(
    refset = refset,
    dcmethod = dcmethod,
    cell_types = cell_types,
    marker_genes = marker_genes,
    phenotype_info = phenotype_info,
    experiment = experiment
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: Single Cell AI Summary
# -----------------------------------------------------------------------------

#' Single Cell AI Summary UI
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

#' Single Cell AI Summary Server
#'
#' Server logic for single-cell AI-generated summaries.
#' Uses omicsai_text_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param refset Reactive returning selected reference dataset
#' @param dcmethod Reactive returning selected deconvolution method
#' @param session Shiny session object (required for getUserOption)
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_text_server
singlecell_ai_summary_server <- function(id,
                                         pgx,
                                         refset,
                                         dcmethod,
                                         session,
                                         cache = NULL,
                                         watermark = FALSE) {
  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    ref <- refset()
    method <- dcmethod()
    shiny::req(ref, method)

    singlecell_build_ai_params(
      pgx = pgx,
      refset = ref,
      dcmethod = method,
      ntop = 20
    )
  })

  # Load template from board.singlecell prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.singlecell/prompts/singlecell_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omicsai::omicsai_load_template(board_template_path)
  })

  # Load single-cell methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.singlecell/prompts/SingleCell_methods.md"
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
