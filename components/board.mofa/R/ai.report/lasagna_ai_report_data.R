##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

LASAGNA_PROMPTS_DIR <- file.path(OPG, "components/board.mofa/prompts/lasagna")

lasagna_ai_md_table <- function(dt) {
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return("No data available.")

  numeric_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
  formatters <- lapply(stats::setNames(numeric_cols, numeric_cols), function(x) {
    function(v) omicsai::omicsai_format_num(v, 3)
  })

  paste(omicsai::omicsai_format_mdtable(dt, formatters = formatters), collapse = "\n")
}

lasagna_ai_build_summary_params <- function(res, contrast, pgx, ntop = 12L) {
  ctx <- lasagna_ai_extract_context(res, contrast, pgx, ntop = ntop)
  if (is.null(ctx)) return(NULL)

  layer_summary <- if (length(ctx$layer_counts) > 0) {
    paste(sprintf("%s=%s", names(ctx$layer_counts), as.integer(ctx$layer_counts)), collapse = ", ")
  } else {
    "No layer information"
  }

  list(
    contrast = ctx$contrast,
    experiment = ctx$experiment,
    headline_metrics = paste0(
      "Nodes: ", ctx$network$n_nodes,
      " | Edges: ", ctx$network$n_edges,
      " | Layers: ", ctx$network$n_layers,
      " | Inter-layer edges: ", ctx$network$n_inter_edges,
      " | Intra-layer edges: ", ctx$network$n_intra_edges
    ),
    layer_summary = layer_summary,
    top_nodes_table = lasagna_ai_md_table(ctx$top_nodes),
    top_edges_table = lasagna_ai_md_table(ctx$top_edges),
    reliability_notes = if (length(ctx$caveats) > 0) {
      paste0("- ", ctx$caveats, collapse = "\n")
    } else {
      "- No major caveats detected for the selected contrast."
    }
  )
}

lasagna_ai_build_report_tables <- function(res, contrast, pgx, ntop = 12L) {
  experiment <- pgx$name %||% pgx$description %||% "omics experiment"

  ctx <- lasagna_ai_extract_context(res, contrast, pgx, ntop = ntop)
  if (is.null(ctx)) {
    return(list(
      experiment = experiment, contrast = contrast %||% "N/A",
      network_snapshot = "", layer_participation = "",
      top_nodes_table = "", top_edges_table = "",
      caveats = "No LASAGNA context available for the selected contrast."
    ))
  }

  ## ---- Network snapshot ----
  network_snapshot <- paste(
    paste0("- Nodes: ", ctx$network$n_nodes),
    paste0("- Edges: ", ctx$network$n_edges),
    paste0("- Layers: ", ctx$network$n_layers),
    paste0("- Inter-layer edges: ", ctx$network$n_inter_edges),
    paste0("- Intra-layer edges: ", ctx$network$n_intra_edges),
    sep = "\n"
  )

  ## ---- Layer participation ----
  layer_participation <- if (length(ctx$layer_counts) > 0) {
    paste(sprintf("- %s: %s", names(ctx$layer_counts), as.integer(ctx$layer_counts)), collapse = "\n")
  } else {
    "No layer labels found in graph."
  }

  ## ---- Caveats ----
  caveats_text <- if (length(ctx$caveats) > 0) {
    paste0("Caveats: ", paste(ctx$caveats, collapse = " "))
  } else {
    "Caveats: none flagged."
  }

  list(
    experiment         = ctx$experiment,
    contrast           = ctx$contrast,
    network_snapshot   = network_snapshot,
    layer_participation = layer_participation,
    top_nodes_table    = lasagna_ai_md_table(ctx$top_nodes),
    top_edges_table    = lasagna_ai_md_table(ctx$top_edges),
    caveats            = caveats_text
  )
}

lasagna_ai_build_methods <- function(pgx, contrast) {
  template <- omicsai::omicsai_load_template(
    file.path(LASAGNA_PROMPTS_DIR, "lasagna_methods.md")
  )

  report_date <- format(Sys.Date(), "%Y-%m-%d")
  params <- list(
    experiment = pgx$name %||% pgx$description %||% "omics experiment",
    contrast = contrast %||% "N/A",
    date = report_date
  )

  omicsai::collapse_lines(
    omicsai::omicsai_substitute_template(template, params),
    sprintf("_This report was generated with OmicsPlayground (BigOmics, %s)._", report_date),
    "_Note: AI-generated interpretation may contain inaccuracies and must be independently verified._",
    sep = "\n\n"
  )
}
