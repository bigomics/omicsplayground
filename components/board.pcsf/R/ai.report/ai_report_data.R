##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

PCSF_PROMPTS_DIR <- file.path(OPG, "components/board.pcsf/prompts")

pcsf_percent <- function(x, digits = 1L) {
  sprintf(paste0("%.", digits, "f%%"), 100 * x)
}

pcsf_compact_hub_table <- function(dt, n = 12L) {
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return("No hub genes available.")
  dt <- head(dt, as.integer(n))

  symbol <- if ("symbol" %in% colnames(dt)) dt$symbol else rownames(dt)
  description <- if ("gene_title" %in% colnames(dt)) dt$gene_title else rep("", nrow(dt))
  node_type <- if ("node_type" %in% colnames(dt)) dt$node_type else rep("", nrow(dt))

  fmt <- data.frame(
    Gene = symbol,
    logFC = dt$logFC,
    Centrality = dt$centrality,
    Type = ifelse(is.na(node_type) | node_type == "", "-", node_type),
    Description = ifelse(is.na(description) | description == "", "-", substr(description, 1, 60)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  omicsai::omicsai_format_mdtable(
    fmt,
    formatters = list(
      logFC = function(x) omicsai::omicsai_format_num(x, 3),
      Centrality = function(x) omicsai::omicsai_format_num(x, 4)
    )
  )
}

pcsf_compact_pathway_table <- function(dt, n = 10L) {
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return("No pathway overlap signal available.")
  dt <- head(dt, as.integer(n))

  fmt <- data.frame(
    Pathway = dt$pathway,
    logFC = dt$logFC,
    `Network genes` = dt$network_genes,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  omicsai::omicsai_format_mdtable(
    fmt,
    formatters = list(
      logFC = function(x) omicsai::omicsai_format_num(x, 3)
    )
  )
}

pcsf_build_summary_params <- function(pgx,
                                      contrast,
                                      pcsf_graph = NULL,
                                      centrality_table = NULL,
                                      ntop = 12L,
                                      pcsf_params = NULL) {
  ctx <- extract_pcsf_context_data(
    pgx = pgx,
    contrast = contrast,
    pcsf_graph = pcsf_graph,
    centrality_table = centrality_table,
    n_top = ntop,
    pcsf_params = pcsf_params
  )
  if (is.null(ctx)) return(NULL)

  n <- ctx$network
  rel <- ctx$reliability

  headline_metrics <- paste0(
    "Nodes: ", n$n_nodes,
    " | edges: ", n$n_edges,
    " | terminal: ", n$n_terminal,
    " | steiner: ", n$n_steiner,
    " | up/down: ", n$up_nodes, "/", n$down_nodes
  )

  reliability_notes <- c(
    paste0("Steiner node fraction: ", pcsf_percent(rel$frac_steiner)),
    paste0("Reportable hub genes: ", rel$n_hubs_reportable),
    paste0("Pathway signal detected: ", ifelse(rel$has_pathway_signal, "yes", "no"))
  )
  if (length(rel$caveats) > 0) reliability_notes <- c(reliability_notes, rel$caveats)

  list(
    contrast = ctx$contrast,
    phenotype = gsub("_", " ", ctx$contrast),
    experiment = ctx$experiment,
    headline_metrics = headline_metrics,
    network_summary = paste(
      paste0("- Total nodes: ", n$n_nodes),
      paste0("- Total edges: ", n$n_edges),
      paste0("- Components: ", n$n_components),
      paste0("- Terminal nodes: ", n$n_terminal),
      paste0("- Steiner nodes: ", n$n_steiner),
      paste0("- Network density: ", omicsai::omicsai_format_num(n$density, 4)),
      sep = "\n"
    ),
    hub_genes_table = pcsf_compact_hub_table(ctx$hubs, n = ntop),
    network_pathways_table = pcsf_compact_pathway_table(ctx$pathways, n = 10L),
    reliability_notes = paste0("- ", reliability_notes, collapse = "\n")
  )
}

pcsf_rank_contexts <- function(context_data) {
  if (length(context_data) == 0) return(data.frame())

  rows <- lapply(context_data, function(ctx) {
    n <- ctx$network
    rel <- ctx$reliability
    max_cent <- if (is.data.frame(ctx$hubs) && nrow(ctx$hubs) > 0 && "centrality" %in% colnames(ctx$hubs)) {
      max(ctx$hubs$centrality, na.rm = TRUE)
    } else {
      0
    }

    score <-
      0.35 * min(n$n_nodes / 120, 1) +
      0.25 * min((n$n_edges + 1) / (n$n_nodes + 1), 1) +
      0.25 * min(rel$n_hubs_reportable / 15, 1) +
      0.15 * as.numeric(rel$has_pathway_signal)

    tier <- if (n$n_nodes < 20 || rel$n_hubs_reportable < 5) {
      "data-limited"
    } else if (score >= 0.70) {
      "strong"
    } else if (score >= 0.45) {
      "moderate"
    } else {
      "weak"
    }

    data.frame(
      contrast = ctx$contrast,
      score = score,
      tier = tier,
      n_nodes = n$n_nodes,
      n_hubs = rel$n_hubs_reportable,
      max_centrality = max_cent,
      has_pathway_signal = rel$has_pathway_signal,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out[order(-out$score), , drop = FALSE]
}

pcsf_build_report_tables <- function(pgx,
                                     contrasts = NULL,
                                     max_contexts = 6L,
                                     ntop = 10L,
                                     ntop_network = 750L,
                                     pcsf_params = NULL) {
  params <- modifyList(
    list(ntop = ntop_network, as_prize = "fc", add_vhce = TRUE, meta_solution = FALSE),
    pcsf_params %||% list()
  )

  all_contrasts <- contrasts %||% pcsf_get_contrasts(pgx)
  if (length(all_contrasts) == 0) {
    return(list(text = "No contrasts available.", data = list(), ranking = data.frame()))
  }

  contexts <- list()
  for (ct in all_contrasts) {
    graph <- pcsf_compute_graph(
      pgx = pgx,
      contrast = ct,
      ntop = params$ntop %||% ntop_network,
      as_prize = params$as_prize %||% "fc",
      add_vhce = isTRUE(params$add_vhce),
      meta_solution = isTRUE(params$meta_solution)
    )
    ctx <- extract_pcsf_context_data(
      pgx = pgx,
      contrast = ct,
      pcsf_graph = graph,
      centrality_table = NULL,
      n_top = ntop,
      pcsf_params = params
    )
    if (!is.null(ctx)) contexts[[ct]] <- ctx
  }

  if (length(contexts) == 0) {
    return(list(text = "No reportable PCSF contexts.", data = list(), ranking = data.frame()))
  }

  ranking <- pcsf_rank_contexts(contexts)
  keep <- head(ranking$contrast, as.integer(max_contexts))
  contexts <- contexts[keep]

  rank_table <- omicsai::omicsai_format_mdtable(
    head(ranking, as.integer(max_contexts)),
    formatters = list(
      score = function(x) omicsai::omicsai_format_num(x, 2),
      max_centrality = function(x) omicsai::omicsai_format_num(x, 4)
    )
  )

  lines <- c(
    paste0("EXPERIMENT: ", pgx$name %||% pgx$description %||% "omics experiment"),
    "BOARD: PCSF",
    "",
    "## Contrast Ranking",
    rank_table,
    "",
    "## Contrast Detail"
  )

  for (ct in keep) {
    ctx <- contexts[[ct]]
    n <- ctx$network
    rel <- ctx$reliability

    lines <- c(
      lines,
      paste0("### ", ct),
      paste0("- Nodes: ", n$n_nodes,
             " | edges: ", n$n_edges,
             " | terminal: ", n$n_terminal,
             " | steiner: ", n$n_steiner,
             " | components: ", n$n_components),
      paste0("- Up/down nodes: ", n$up_nodes, "/", n$down_nodes,
             " | pathway signal: ", ifelse(rel$has_pathway_signal, "yes", "no"),
             " | steiner fraction: ", pcsf_percent(rel$frac_steiner)),
      "",
      "Top hub genes:",
      pcsf_compact_hub_table(ctx$hubs, n = ntop),
      "",
      "Pathway overlap:",
      pcsf_compact_pathway_table(ctx$pathways, n = ntop),
      "",
      if (length(rel$caveats) > 0) paste0("Caveats: ", paste(rel$caveats, collapse = " ")) else "Caveats: none flagged.",
      ""
    )
  }

  list(
    text = paste(lines, collapse = "\n"),
    data = contexts,
    ranking = ranking
  )
}

pcsf_build_methods <- function(pgx) {
  template <- omicsai::omicsai_load_template(
    file.path(PCSF_PROMPTS_DIR, "pcsf_report_methods.md")
  )

  params <- list(
    experiment = pgx$name %||% pgx$description %||% "omics experiment",
    date = as.character(Sys.Date())
  )

  omicsai::omicsai_substitute_template(template, params)
}
