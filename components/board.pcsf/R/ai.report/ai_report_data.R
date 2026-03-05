##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

PCSF_PROMPTS_DIR <- file.path(OPG, "components/board.pcsf/prompts")

pcsf_percent <- function(x, digits = 1L) {
  sprintf(paste0("%.", digits, "f%%"), 100 * x)
}

pcsf_yes_no <- function(x) {
  if (isTRUE(x)) "yes" else "no"
}

pcsf_format_num_safe <- function(x, digits = 3L) {
  if (length(x) == 0 || is.na(x)) return("NA")
  omicsai::omicsai_format_num(x, digits)
}

pcsf_format_qv <- function(x) {
  if (length(x) == 0 || is.na(x)) return("NA")
  formatC(x, format = "e", digits = 1)
}

pcsf_detect_order_column <- function(samples) {
  if (is.null(samples) || !is.data.frame(samples) || ncol(samples) == 0) return(NULL)
  cols <- colnames(samples)
  lower_cols <- tolower(cols)
  match_idx <- match(c("time", "timepoint", "day", "hour", "stage"), lower_cols, nomatch = 0L)
  match_idx <- match_idx[match_idx > 0L]
  if (length(match_idx) == 0) return(NULL)
  cols[match_idx[[1]]]
}

pcsf_build_dataset_context_md <- function(pgx, contrasts) {
  sample_cols <- if (!is.null(pgx$samples) && is.data.frame(pgx$samples)) colnames(pgx$samples) else character(0)
  order_col <- pcsf_detect_order_column(pgx$samples)

  lines <- c(
    "## Dataset Context",
    paste0("- Experiment: ", pgx$name %||% "omics experiment"),
    paste0("- Organism: ", pgx$organism %||% "unknown"),
    paste0("- Data type: ", pgx$datatype %||% "unknown"),
    paste0("- Included contrasts: ", paste(contrasts, collapse = ", "))
  )

  if (!is.null(pgx$description) && nzchar(pgx$description) && !identical(pgx$description, pgx$name)) {
    lines <- c(lines, paste0("- Dataset description: ", pgx$description))
  }
  if (length(sample_cols) > 0) {
    lines <- c(lines, paste0("- Available sample annotation columns: ", paste(sample_cols, collapse = ", ")))
  }
  if (!is.null(order_col)) {
    lines <- c(lines, paste0("- Ordered-condition hint from sample annotations: ", order_col))
  }

  lines
}

pcsf_build_fact_constraints_md <- function(contexts) {
  lines <- c("## Per-Contrast Factual Constraints")

  for (ct in names(contexts)) {
    ctx <- contexts[[ct]]
    rel <- ctx$reliability
    lines <- c(
      lines,
      paste0(
        "- ", ct,
        ": steiner nodes = ", ctx$network$n_steiner,
        "; steiner fraction = ", pcsf_percent(rel$frac_steiner),
        "; pathway signal = ", pcsf_yes_no(rel$has_pathway_signal),
        "; top-hub effects all |logFC| < 0.5 = ", pcsf_yes_no(rel$hubs_have_modest_effects)
      )
    )
  }

  lines
}

pcsf_build_reference_catalog <- function(contexts,
                                         max_refs = 15L,
                                         top_pathways_per_contrast = 1L,
                                         max_recurrent_pathways = 4L) {
  rows <- list()
  idx <- 1L

  add_row <- function(text, group, priority) {
    rows[[length(rows) + 1L]] <<- data.frame(
      text = text,
      group = as.integer(group),
      priority = as.numeric(priority),
      order = idx,
      stringsAsFactors = FALSE
    )
    idx <<- idx + 1L
  }

  for (i in seq_along(contexts)) {
    ct <- names(contexts)[[i]]
    ctx <- contexts[[i]]
    rel <- ctx$reliability

    metric_bits <- c(
      paste0("nodes = ", ctx$network$n_nodes),
      paste0("steiner = ", ctx$network$n_steiner),
      paste0("steiner fraction = ", pcsf_percent(rel$frac_steiner)),
      paste0("pathway signal = ", pcsf_yes_no(rel$has_pathway_signal))
    )
    if (isTRUE(rel$hubs_have_modest_effects)) {
      metric_bits <- c(metric_bits, "top hub |logFC| values all < 0.5")
    }

    add_row(
      paste0(ct, " — network metrics: ", paste(metric_bits, collapse = "; ")),
      group = 1L,
      priority = i
    )
  }

  for (i in seq_along(contexts)) {
    ct <- names(contexts)[[i]]
    pathways <- contexts[[i]]$pathways
    if (is.null(pathways) || !is.data.frame(pathways) || nrow(pathways) == 0) next

    pathways <- pathways[order(-abs(pathways$logFC), -pathways$network_genes, pathways$pathway), , drop = FALSE]
    pathways <- head(pathways, as.integer(top_pathways_per_contrast))
    for (j in seq_len(nrow(pathways))) {
      p <- pathways[j, ]
      qv_txt <- if ("qv" %in% colnames(pathways) && !is.na(p$qv)) {
        paste0("; q = ", pcsf_format_qv(p$qv))
      } else {
        ""
      }
      add_row(
        paste0(
          ct, " — ", p$pathway,
          " (pathway logFC = ", pcsf_format_num_safe(p$logFC, 3L),
          "; overlap = ", p$network_genes,
          qv_txt, ")"
        ),
        group = 2L,
        priority = i + (j / 100)
      )
    }
  }

  pathway_rows <- list()
  for (ct in names(contexts)) {
    ctx <- contexts[[ct]]
    if (is.null(ctx$pathways) || !is.data.frame(ctx$pathways) || nrow(ctx$pathways) == 0) next
    dt <- ctx$pathways
    dt$contrast <- ct
    pathway_rows[[length(pathway_rows) + 1L]] <- dt
  }

  if (length(pathway_rows) > 0) {
    pathway_dt <- do.call(rbind, pathway_rows)
    split_idx <- split(seq_len(nrow(pathway_dt)), pathway_dt$pathway)
    recurrent <- Filter(function(ix) length(unique(pathway_dt$contrast[ix])) >= 2, split_idx)

    if (length(recurrent) > 0) {
      recurrent_rows <- lapply(names(recurrent), function(pathway_name) {
        ix <- recurrent[[pathway_name]]
        subset <- pathway_dt[ix, , drop = FALSE]
        data.frame(
          pathway = pathway_name,
          n_contrasts = length(unique(subset$contrast)),
          max_abs_logFC = max(abs(subset$logFC), na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      })
      recurrent_rows <- do.call(rbind, recurrent_rows)
      recurrent_rows <- recurrent_rows[order(-recurrent_rows$n_contrasts, -recurrent_rows$max_abs_logFC, recurrent_rows$pathway), , drop = FALSE]
      recurrent_rows <- head(recurrent_rows, as.integer(max_recurrent_pathways))

      for (i in seq_len(nrow(recurrent_rows))) {
        subset <- pathway_dt[pathway_dt$pathway == recurrent_rows$pathway[i], , drop = FALSE]
        subset <- subset[order(subset$contrast), , drop = FALSE]
        qv_vals <- subset$qv[!is.na(subset$qv)]
        qv_txt <- if (length(qv_vals) > 0) {
          paste0("; q range = ", pcsf_format_qv(min(qv_vals, na.rm = TRUE)), " to ", pcsf_format_qv(max(qv_vals, na.rm = TRUE)))
        } else {
          ""
        }
        add_row(
          paste0(
            recurrent_rows$pathway[i],
            " — recurring pathway across ",
            recurrent_rows$n_contrasts[i],
            " contrasts (",
            paste(subset$contrast, collapse = ", "),
            "); pathway logFC = ",
            pcsf_format_num_safe(min(subset$logFC, na.rm = TRUE), 3L),
            " to ",
            pcsf_format_num_safe(max(subset$logFC, na.rm = TRUE), 3L),
            "; overlap = ",
            min(subset$network_genes, na.rm = TRUE),
            " to ",
            max(subset$network_genes, na.rm = TRUE),
            qv_txt
          ),
          group = 3L,
          priority = i
        )
      }
    }
  }

  if (length(rows) == 0) return(data.frame())

  ref_df <- do.call(rbind, rows)
  ref_df <- ref_df[order(ref_df$group, ref_df$priority, ref_df$order), , drop = FALSE]
  ref_df <- head(ref_df, as.integer(max_refs))
  ref_df$ref_no <- seq_len(nrow(ref_df))
  ref_df
}

pcsf_render_reference_catalog_md <- function(ref_df) {
  lines <- c(
    "## Referenceable Evidence",
    "Use only these pre-numbered [n] entries when citing evidence in the narrative.",
    "The final `## Data References` section should copy only the cited entries below, preserving [n] numbers and wording exactly."
  )

  if (is.null(ref_df) || !is.data.frame(ref_df) || nrow(ref_df) == 0) {
    return(c(lines, "- No referenceable evidence catalog available."))
  }

  c(lines, paste0("- [", ref_df$ref_no, "] ", ref_df$text))
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

  paste(omicsai::omicsai_format_mdtable(
    fmt,
    formatters = list(
      logFC = function(x) omicsai::omicsai_format_num(x, 3),
      Centrality = function(x) omicsai::omicsai_format_num(x, 4)
    )
  ), collapse = "\n")
}

pcsf_compact_pathway_table <- function(dt, n = 10L) {
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return("No pathway overlap signal available.")
  dt <- head(dt, as.integer(n))

  fmt <- data.frame(
    Pathway = dt$pathway,
    logFC = dt$logFC,
    q = dt$qv,
    `Network genes` = dt$network_genes,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  paste(omicsai::omicsai_format_mdtable(
    fmt,
    formatters = list(
      logFC = function(x) omicsai::omicsai_format_num(x, 3),
      q = function(x) ifelse(is.na(x), "NA", pcsf_format_qv(x))
    )
  ), collapse = "\n")
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
    paste0("Pathway signal detected: ", ifelse(rel$has_pathway_signal, "yes", "no")),
    paste0("Top hub effects all |logFC| < 0.5: ", pcsf_yes_no(rel$hubs_have_modest_effects))
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
  ref_catalog <- pcsf_build_reference_catalog(contexts, max_refs = 15L)

  rank_table <- omicsai::omicsai_format_mdtable(
    head(ranking, as.integer(max_contexts)),
    formatters = list(
      score = function(x) omicsai::omicsai_format_num(x, 2),
      max_centrality = function(x) omicsai::omicsai_format_num(x, 4)
    )
  )

  lines <- c(
    pcsf_build_dataset_context_md(pgx, keep),
    "",
    pcsf_build_fact_constraints_md(contexts),
    "",
    pcsf_render_reference_catalog_md(ref_catalog),
    "",
    "## Raw Report Data",
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
    file.path(PCSF_PROMPTS_DIR, "pcsf_methods.md")
  )

  params <- list(
    experiment = pgx$name %||% pgx$description %||% "omics experiment",
    date = as.character(Sys.Date())
  )

  omicsai::omicsai_substitute_template(template, params)
}
