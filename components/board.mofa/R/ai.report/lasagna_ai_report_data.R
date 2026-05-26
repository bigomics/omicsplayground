##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# LASAGNA AI Report — data extraction and template rendering
# =============================================================================
# All deterministic data functions consumed by the AI-report Shiny module
# (lasagna_ai_text_server.R). No Shiny dependencies. No LLM calls.
#
# Layout (top to bottom):
#   - verbalisers (board-local; consolidated for the future compaction pass)
#   - module decomposition helpers (.lasagna_compute_modules, .lasagna_tier)
#   - leaf renderers, ONE per {{placeholder}} in lasagna_report_data.md:
#       lasagna_data_overview, lasagna_data_contrast_block,
#       lasagna_data_layer_participation, lasagna_data_modules_summary,
#       lasagna_data_module_bridges, lasagna_data_module_detail
#   - per-module block renderer:  .lasagna_render_module_block (MODULE_TEMPLATE)
#   - orchestrators:              lasagna_build_report_tables,
#                                 lasagna_build_summary_params,
#                                 lasagna_build_methods
#
# Authoritative thresholds also appear in lasagna_interpretation.md (any
# change there MUST be mirrored here and vice versa — the LLM is told
# what each label means).

LASAGNA_PROMPTS_DIR <- file.path(OPG, "components/board.mofa/prompts/lasagna")


# =============================================================================
# Verbalisers — keep together; promote to omicsai during the compaction pass
# =============================================================================

#' Verbalise a cross-view network centrality score in [0, 1].
#' Bins: hub (≥ 0.8) / central (≥ 0.6) / intermediate (≥ 0.3) / peripheral.
lasagna_verbalize_centrality <- function(c_, breaks = c(0.3, 0.6, 0.8),
                                         na_label = "not tested") {
  out <- rep(na_label, length(c_))
  ok  <- !is.na(c_)
  out[ok] <- ifelse(c_[ok] >= breaks[3], "hub",
              ifelse(c_[ok] >= breaks[2], "central",
              ifelse(c_[ok] >= breaks[1], "intermediate", "peripheral")))
  out
}

#' Verbalise a per-layer participation fraction (0–1) within a module.
#' Bins: dominant (≥ 60%) / major (≥ 30%) / minor (≥ 10%) / marginal.
lasagna_verbalize_layer_pct <- function(p, breaks = c(0.10, 0.30, 0.60),
                                        na_label = "not tested") {
  out <- rep(na_label, length(p))
  ok  <- !is.na(p)
  out[ok] <- ifelse(p[ok] >= breaks[3], "dominant",
              ifelse(p[ok] >= breaks[2], "major",
              ifelse(p[ok] >= breaks[1], "minor", "marginal")))
  out
}

#' Verbalise the cross-layer connectivity of a module.
#' Bins: rich-bridged (≥ 10 inter-layer edges) / well-bridged (≥ 5) /
#'       lightly-bridged (≥ 1) / single-layer.
lasagna_verbalize_connectivity <- function(n_inter,
                                           breaks = c(1, 5, 10),
                                           na_label = "not tested") {
  out <- rep(na_label, length(n_inter))
  ok  <- !is.na(n_inter)
  out[ok] <- ifelse(n_inter[ok] >= breaks[3], "rich-bridged",
              ifelse(n_inter[ok] >= breaks[2], "well-bridged",
              ifelse(n_inter[ok] >= breaks[1], "lightly-bridged",
                                                "single-layer")))
  out
}

#' Look up one-line functional descriptions for features. NULL-safe.
#' Local to LASAGNA (board-prefixed) to avoid the collision footgun
#' that bit us on data_overview — promote to a shared helper during
#' the compaction pass.
.lasagna_resolve_functions <- function(features, annot, max_chars = 60L) {
  funcs <- rep("", length(features))
  if (is.null(annot) || length(features) == 0) return(funcs)
  func_col <- intersect(c("gene_title", "gene_name", "description"),
                        colnames(annot))
  if (length(func_col) == 0) return(funcs)
  idx <- match(features, rownames(annot))
  valid <- !is.na(idx)
  funcs[valid] <- as.character(annot[idx[valid], func_col[1]])
  substr(funcs, 1, max_chars)
}

#' Verbalise whole-graph density.
#' Bins: dense (≥ 0.1) / moderate (≥ 0.02) / sparse (≥ 0.005) / fragmented.
lasagna_verbalize_density <- function(d, breaks = c(0.005, 0.02, 0.10),
                                      na_label = "not tested") {
  out <- rep(na_label, length(d))
  ok  <- !is.na(d)
  out[ok] <- ifelse(d[ok] >= breaks[3], "dense",
              ifelse(d[ok] >= breaks[2], "moderate",
              ifelse(d[ok] >= breaks[1], "sparse", "fragmented")))
  out
}


# =============================================================================
# Module decomposition (community detection on the multipartite graph)
# =============================================================================

#' Decompose the LASAGNA graph into modules via Louvain community detection.
#'
#' Returns a list keyed `M1`, `M2`, … (descending by tier score). Each
#' element carries: node ids, per-layer counts, intra/inter edge counts,
#' top hub nodes (by centrality), top inter-layer edges (by weight).
.lasagna_compute_modules <- function(graph, nodes, edges,
                                     ntop_nodes = 10L,
                                     ntop_edges = 10L,
                                     min_size = 3L) {
  if (is.null(graph) || !inherits(graph, "igraph")) {
    message("[lasagna] compute_modules: graph is NULL or not igraph")
    return(list())
  }
  if (igraph::vcount(graph) < min_size) {
    message("[lasagna] compute_modules: vcount(", igraph::vcount(graph),
            ") < min_size(", min_size, ")")
    return(list())
  }

  undirected <- igraph::as.undirected(graph, mode = "collapse")

  ## Louvain requires non-negative edge weights. LASAGNA edges carry
  ## signed `rho` (and a `weight` attribute that may inherit the sign),
  ## so the default weighted Louvain silently errors and we lose all
  ## community structure. Use absolute weights when present, NA when
  ## absent (NA = unweighted Louvain).
  lw <- igraph::edge_attr(undirected, "weight")
  louvain_weights <- if (!is.null(lw)) abs(as.numeric(lw)) else NA
  comm <- tryCatch(igraph::cluster_louvain(undirected, weights = louvain_weights),
                   error = function(e) {
                     message("[lasagna] cluster_louvain failed: ",
                             conditionMessage(e))
                     NULL
                   })
  if (is.null(comm)) return(list())

  centrality <- tryCatch(
    igraph::eigen_centrality(undirected, scale = TRUE,
                             weights = louvain_weights)$vector,
    error = function(e) {
      message("[lasagna] eigen_centrality failed: ", conditionMessage(e))
      setNames(rep(0, igraph::vcount(graph)), igraph::V(graph)$name)
    }
  )

  membership <- igraph::membership(comm)
  size_tab <- table(membership)
  keep_ids <- as.integer(names(size_tab)[as.integer(size_tab) >= min_size])
  if (length(keep_ids) == 0) {
    message("[lasagna] compute_modules: Louvain found ", length(size_tab),
            " communities but none reach min_size(", min_size,
            "); top sizes = ",
            paste(head(sort(as.integer(size_tab), decreasing = TRUE), 6),
                  collapse = ","))
    return(list())
  }

  node_layer <- setNames(as.character(nodes$layer),
                         as.character(nodes$id))
  from_layer <- node_layer[as.character(edges$from)]
  to_layer   <- node_layer[as.character(edges$to)]
  edge_inter <- !is.na(from_layer) & !is.na(to_layer) &
                from_layer != to_layer

  out <- list()
  for (cid in keep_ids) {
    member_ids <- names(membership)[membership == cid]
    layer_counts <- sort(table(node_layer[member_ids]), decreasing = TRUE)
    n_layers_spanned <- length(layer_counts)

    e_idx <- which(as.character(edges$from) %in% member_ids &
                   as.character(edges$to)   %in% member_ids)
    n_intra <- sum(!edge_inter[e_idx])
    n_inter <- sum(edge_inter[e_idx])

    cent <- centrality[member_ids]
    hub_order <- order(-cent, na.last = TRUE)
    top_n <- head(member_ids[hub_order], ntop_nodes)
    node_idx <- match(top_n, as.character(nodes$id))
    top_nodes <- data.frame(
      id         = top_n,
      label      = as.character(nodes$label[node_idx]),
      layer      = node_layer[top_n],
      centrality = as.numeric(cent[top_n]),
      rho        = if ("rho" %in% colnames(nodes)) as.numeric(nodes$rho[node_idx])
                   else rep(NA_real_, length(top_n)),
      fc         = if ("fc"  %in% colnames(nodes)) as.numeric(nodes$fc[node_idx])
                   else rep(NA_real_, length(top_n)),
      stringsAsFactors = FALSE
    )

    ## Module-level trait association: aggregate rho across ALL module
    ## members (not just top hubs) so a strong signal in non-hub
    ## members is not lost. Signed by the majority direction.
    member_idx <- match(member_ids, as.character(nodes$id))
    mod_rho <- if ("rho" %in% colnames(nodes)) {
      as.numeric(nodes$rho[member_idx])
    } else NA_real_
    rho_summary <- if (any(!is.na(mod_rho))) {
      sgn <- sign(sum(mod_rho, na.rm = TRUE))
      if (is.na(sgn) || sgn == 0) sgn <- 1
      mean(abs(mod_rho), na.rm = TRUE) * sgn
    } else NA_real_

    e_member_inter <- e_idx[edge_inter[e_idx]]
    if (length(e_member_inter) > 0) {
      w <- if ("rho" %in% colnames(edges)) edges$rho[e_member_inter]
           else if ("weight" %in% colnames(edges)) edges$weight[e_member_inter]
           else rep(NA_real_, length(e_member_inter))
      ord <- order(-abs(w), na.last = TRUE)
      pick <- head(e_member_inter[ord], ntop_edges)
      top_edges <- data.frame(
        from_layer = from_layer[pick],
        from_node  = as.character(edges$from[pick]),
        to_layer   = to_layer[pick],
        to_node    = as.character(edges$to[pick]),
        rho        = if ("rho" %in% colnames(edges)) edges$rho[pick]
                     else rep(NA_real_, length(pick)),
        stringsAsFactors = FALSE
      )
    } else {
      top_edges <- data.frame(
        from_layer = character(0), from_node = character(0),
        to_layer   = character(0), to_node   = character(0),
        rho = numeric(0), stringsAsFactors = FALSE
      )
    }

    ## Cross-layer hubs: top hubs that touch at least one inter-layer edge.
    inter_endpoints <- unique(c(
      as.character(edges$from[e_idx[edge_inter[e_idx]]]),
      as.character(edges$to[  e_idx[edge_inter[e_idx]]])
    ))
    cross_layer_hubs <- intersect(top_n, inter_endpoints)

    out[[as.character(cid)]] <- list(
      members          = member_ids,
      n_nodes          = length(member_ids),
      layer_counts     = layer_counts,
      n_layers_spanned = n_layers_spanned,
      n_intra_edges    = n_intra,
      n_inter_edges    = n_inter,
      top_nodes        = top_nodes,
      top_edges        = top_edges,
      cross_layer_hubs = cross_layer_hubs,
      rho_summary      = rho_summary
    )
  }

  ## Tier score: rewards multi-layer membership + inter-layer bridges +
  ## node count. Used to order modules descending.
  scores <- vapply(out, function(m) {
    0.45 * min(m$n_inter_edges / 10, 1) +
    0.30 * min(m$n_layers_spanned / 3, 1) +
    0.25 * min(m$n_nodes / 30, 1)
  }, numeric(1))
  out <- out[order(-scores)]
  names(out) <- paste0("M", seq_along(out))
  message("[lasagna] compute_modules: returned ", length(out),
          " modules with sizes = ",
          paste(vapply(out, function(m) m$n_nodes, integer(1)),
                collapse = ","))
  out
}

#' For one focal module, return its top partner modules by shared
#' membership and shared inter-layer edges. Mirrors the multiwgcna
#' "module ↔ partner-module" cross-layer linkage pattern but uses
#' graph membership (LASAGNA has no eigengenes).
#'
#' Score = (shared hub count) · 2 + (shared member count) / 10. The
#' weight on shared hubs keeps the LLM focused on biologically
#' meaningful overlap rather than incidental membership.
.lasagna_module_partners <- function(modules, focal_name, top_n = 3L) {
  if (length(modules) < 2L) return(data.frame())
  focal <- modules[[focal_name]]
  others <- setdiff(names(modules), focal_name)

  rows <- list()
  focal_hubs <- if (!is.null(focal$top_nodes)) focal$top_nodes$id else character(0)
  for (om_name in others) {
    om <- modules[[om_name]]
    om_hubs <- if (!is.null(om$top_nodes)) om$top_nodes$id else character(0)
    shared_hubs <- intersect(focal_hubs, om_hubs)
    shared_members <- intersect(focal$members, om$members)
    score <- length(shared_hubs) * 2 + length(shared_members) / 10
    if (score <= 0) next
    rows[[length(rows) + 1L]] <- data.frame(
      partner        = om_name,
      shared_hubs    = length(shared_hubs),
      shared_members = length(shared_members),
      score          = score,
      shared_labels  = paste(shared_hubs, collapse = ","),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(data.frame())
  df <- do.call(rbind, rows)
  df <- df[order(-df$score), , drop = FALSE]
  head(df, as.integer(top_n))
}

#' Tier per module (strong/moderate/weak) from per-module metrics.
.lasagna_tier <- function(m) {
  score <- 0.45 * min(m$n_inter_edges / 10, 1) +
           0.30 * min(m$n_layers_spanned / 3, 1) +
           0.25 * min(m$n_nodes / 30, 1)
  if (score >= 0.60) "strong"
  else if (score >= 0.35) "moderate"
  else "weak"
}


# =============================================================================
# Leaf renderers — ONE per {{placeholder}} in lasagna_report_data.md
# =============================================================================

#' {{experiment}}, {{organism}}, {{contrast}}, {{n_samples}}, {{n_layers}},
#' {{layer_names}}, {{n_nodes}}, {{n_edges}}, {{inter_layer_pct}},
#' {{n_modules_total}}, {{n_modules_used}} — the Overview metadata block.
lasagna_data_overview <- function(ctx, modules, n_modules_used) {
  pgx <- ctx$pgx
  experiment <- multiomics_ai_experiment_label(pgx)
  organism   <- pgx$organism %||% "unknown"
  n_samples  <- if (!is.null(pgx$samples)) nrow(pgx$samples) else NA_integer_

  layer_names <- character(0)
  if (!is.null(ctx$layer_counts)) {
    layer_names <- names(ctx$layer_counts)
  }
  n_layers <- length(layer_names)

  inter_layer_pct <- if (!is.null(ctx$network$n_edges) &&
                         ctx$network$n_edges > 0) {
    sprintf("%s%% (%d / %d)",
            omicsai::omicsai_format_num(
              100 * ctx$network$n_inter_edges / ctx$network$n_edges, 1),
            ctx$network$n_inter_edges, ctx$network$n_edges)
  } else "not available"

  list(
    experiment       = experiment,
    organism         = organism,
    contrast         = ctx$contrast %||% "N/A",
    n_samples        = as.character(n_samples %||% "unknown"),
    n_layers         = as.character(n_layers),
    layer_names      = if (n_layers > 0) paste(layer_names, collapse = ", ")
                       else "unknown",
    n_nodes          = as.character(ctx$network$n_nodes %||% "0"),
    n_edges          = as.character(ctx$network$n_edges %||% "0"),
    inter_layer_pct  = inter_layer_pct,
    n_modules_total  = as.character(length(modules)),
    n_modules_used   = as.character(n_modules_used)
  )
}

#' {{contrast_block}} — contrast definition with sample-group counts when
#' available.
lasagna_data_contrast_block <- function(pgx, contrast) {
  if (is.null(contrast) || !nzchar(contrast)) return("(no contrast selected)")
  ct_def <- tryCatch(playbase::pgx.getContrasts(pgx, contrast),
                     error = function(e) NULL)
  if (is.null(ct_def)) return(sprintf("- %s", contrast))
  if (is.data.frame(ct_def) || is.matrix(ct_def)) {
    tab <- table(as.character(ct_def[, 1]))
    lines <- sprintf("- %s: %s (n=%d)", contrast, names(tab),
                     as.integer(tab))
    return(paste(lines, collapse = "\n"))
  }
  sprintf("- %s", contrast)
}

#' {{layer_participation}} — per-layer node count (raw), verbalised
#' participation share across the whole network.
lasagna_data_layer_participation <- function(ctx) {
  lc <- ctx$layer_counts
  if (is.null(lc) || length(lc) == 0) return("(no layer information)")
  total <- sum(as.integer(lc))
  pct   <- as.integer(lc) / total
  labels <- lasagna_verbalize_layer_pct(pct)
  paste(sprintf("- %s: %d nodes (%s)", names(lc),
                as.integer(lc), labels),
        collapse = "\n")
}

#' {{modules_summary_table}}, {{lead_module}} — cross-module summary
#' table + lead-module identifier (template owns the prose framing).
lasagna_data_modules_summary <- function(modules) {
  if (length(modules) == 0) {
    return(list(table = "(no modules detected)", lead_module = "—"))
  }
  rows <- lapply(names(modules), function(mn) {
    m <- modules[[mn]]
    top_hub <- if (nrow(m$top_nodes) > 0) {
      paste0("*", m$top_nodes$label[1] %||% m$top_nodes$id[1], "*")
    } else "—"
    data.frame(
      Module                  = mn,
      Tier                    = .lasagna_tier(m),
      Nodes                   = as.character(m$n_nodes),
      `Layers spanned`        = as.character(m$n_layers_spanned),
      `Cross-layer edges`     = as.character(m$n_inter_edges),
      `Top hub`               = top_hub,
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  list(table       = paste(omicsai::omicsai_format_mdtable(df), collapse = "\n"),
       lead_module = names(modules)[1])
}

#' {{module_bridges}} — module pairs sharing one or more hub nodes.
lasagna_data_module_bridges <- function(modules) {
  if (length(modules) < 2) return("(too few modules for bridges)")
  hubs <- lapply(modules, function(m) {
    if (nrow(m$top_nodes) > 0) m$top_nodes$id else character(0)
  })
  pairs <- character(0)
  nm <- names(modules)
  for (i in seq_len(length(modules) - 1L)) {
    for (j in (i + 1L):length(modules)) {
      shared <- intersect(hubs[[i]], hubs[[j]])
      if (length(shared) > 0) {
        shared_labels <- modules[[i]]$top_nodes$label[
          match(shared, modules[[i]]$top_nodes$id)] %||% shared
        pairs <- c(pairs, sprintf("- %s ↔ %s: %d shared hub(s) (%s)",
                                  nm[i], nm[j], length(shared),
                                  paste(sprintf("*%s*", shared_labels),
                                        collapse = ", ")))
      }
    }
  }
  if (length(pairs) == 0) return("(no shared-hub bridges between modules)")
  paste(pairs, collapse = "\n")
}

#' {{module_detail}} — concatenated per-module blocks.
#'
#' `annot` is `pgx$genes` (used to look up the `Function` column for
#' top hubs). `contrast` is the active contrast name (used to phrase
#' the per-module trait-association line). `all_modules` is the
#' complete set of modules in the report — passed so each block can
#' compute its top partner modules.
lasagna_data_module_detail <- function(modules,
                                       annot = NULL,
                                       contrast = NULL,
                                       ntop = 10L) {
  if (length(modules) == 0) return("")
  blocks <- vapply(names(modules), function(mn) {
    .lasagna_render_module_block(modules[[mn]], mn,
                                 all_modules = modules,
                                 annot = annot,
                                 contrast = contrast,
                                 ntop = ntop)
  }, character(1))
  paste(blocks, collapse = "\n\n")
}


# =============================================================================
# Per-module block renderer (shared by Summary + Report)
# =============================================================================

## MODULE_TEMPLATE — loaded from disk. Edit prompts/lasagna/lasagna_module_data.md
## directly to change wording; the R renderer only fills values.
MODULE_TEMPLATE_LASAGNA <- omicsai::omicsai_load_template(
  file.path(LASAGNA_PROMPTS_DIR, "lasagna_module_data.md")
)

#' Render one module's data dict into a MODULE_TEMPLATE block. All
#' numeric quantities (centrality, ρ, layer share) are verbalised via
#' the verbalisers at the top of this file.
#'
#' `all_modules` (optional) lets the block include the focal module's
#' top partner modules — the LASAGNA analogue of the multiwgcna
#' "cross-layer module links" section. `annot` (optional) drives the
#' hub-function annotation column; pass `pgx$genes` from the caller.
.lasagna_render_module_block <- function(m, module_name,
                                         all_modules = NULL,
                                         annot = NULL,
                                         contrast = NULL,
                                         ntop = 10L) {
  tier <- .lasagna_tier(m)

  layer_mix <- if (length(m$layer_counts) > 0) {
    paste(sprintf("%s=%d", names(m$layer_counts),
                  as.integer(m$layer_counts)),
          collapse = ", ")
  } else "unknown"

  connectivity_label <- lasagna_verbalize_connectivity(m$n_inter_edges)

  ## Trait association line (per-module rho aggregate, verbalised).
  trait_summary <- if (!is.null(m$rho_summary) && !is.na(m$rho_summary)) {
    sprintf("%s with %s",
            omicsai::omicsai_verbalize_r(m$rho_summary),
            contrast %||% "contrast")
  } else "—"

  nodes_str <- if (!is.null(m$top_nodes) && nrow(m$top_nodes) > 0) {
    n <- head(m$top_nodes, ntop)
    cross_layer <- n$id %in% m$cross_layer_hubs
    funcs <- .lasagna_resolve_functions(n$id, annot)
    funcs[funcs == ""] <- "—"
    df <- data.frame(
      Symbol         = paste0("*", ifelse(is.na(n$label) | n$label == "",
                                          n$id, n$label), "*"),
      Layer          = n$layer,
      `Network role` = lasagna_verbalize_centrality(n$centrality),
      `Cross-layer`  = ifelse(cross_layer, "yes", "no"),
      Function       = funcs,
      check.names = FALSE, stringsAsFactors = FALSE
    )
    paste(omicsai::omicsai_format_mdtable(df), collapse = "\n")
  } else "—"

  edges_str <- if (!is.null(m$top_edges) && nrow(m$top_edges) > 0) {
    e <- head(m$top_edges, ntop)
    df <- data.frame(
      `From layer`     = e$from_layer,
      `From node`      = paste0("*", e$from_node, "*"),
      `To layer`       = e$to_layer,
      `To node`        = paste0("*", e$to_node, "*"),
      `Edge strength`  = omicsai::omicsai_verbalize_r(e$rho),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    paste(omicsai::omicsai_format_mdtable(df), collapse = "\n")
  } else "—"

  cross_hubs_str <- if (length(m$cross_layer_hubs) > 0) {
    labels <- m$top_nodes$label[match(m$cross_layer_hubs, m$top_nodes$id)]
    labels[is.na(labels) | labels == ""] <-
      m$cross_layer_hubs[is.na(labels) | labels == ""]
    paste(sprintf("*%s*", labels), collapse = ", ")
  } else "—"

  ## Per-module partner modules (cross-module linkage).
  partners_str <- if (!is.null(all_modules) && length(all_modules) > 1L) {
    partners <- .lasagna_module_partners(all_modules, module_name, top_n = 3L)
    if (nrow(partners) > 0) {
      paste(sprintf("- %s (%d shared hub(s), %d shared member(s))",
                    partners$partner,
                    partners$shared_hubs,
                    partners$shared_members),
            collapse = "\n")
    } else "—"
  } else "—"

  omicsai::omicsai_substitute_template(MODULE_TEMPLATE_LASAGNA, list(
    module             = module_name,
    n_nodes            = as.character(m$n_nodes),
    tier               = tier,
    layer_mix          = layer_mix,
    connectivity_label = connectivity_label,
    n_inter_edges      = as.character(m$n_inter_edges),
    n_layers_spanned   = as.character(m$n_layers_spanned),
    n_top              = as.character(min(ntop, m$n_nodes)),
    n_top_edges        = as.character(min(ntop, m$n_inter_edges)),
    trait_summary      = trait_summary,
    top_nodes_table    = nodes_str,
    top_edges_table    = edges_str,
    cross_layer_hubs   = cross_hubs_str,
    partner_modules    = partners_str
  ))
}


# =============================================================================
# Orchestrators
# =============================================================================

#' Build structured report tables from LASAGNA results.
#'
#' @return list(text = character, modules = list, ctx = list)
lasagna_build_report_tables <- function(res, contrast, pgx,
                                        n_modules = 8L, ntop = 10L) {
  ctx <- lasagna_ai_extract_context(res, contrast, pgx, ntop = ntop)
  if (is.null(ctx)) {
    return(list(text = "(no LASAGNA context available)",
                modules = list(), ctx = NULL))
  }
  ctx$pgx <- pgx

  modules <- .lasagna_compute_modules(res$graph, ctx$nodes_all,
                                      ctx$edges_all,
                                      ntop_nodes = ntop,
                                      ntop_edges = ntop)
  keep <- head(modules, as.integer(n_modules))

  overview <- lasagna_data_overview(ctx, modules,
                                    n_modules_used = length(keep))
  contrast_block <- lasagna_data_contrast_block(pgx, contrast)
  layer_part     <- lasagna_data_layer_participation(ctx)
  modsum         <- lasagna_data_modules_summary(keep)
  bridges        <- lasagna_data_module_bridges(keep)
  module_detail  <- lasagna_data_module_detail(keep,
                                               annot = pgx$genes,
                                               contrast = contrast,
                                               ntop = ntop)

  tmpl <- omicsai::omicsai_load_template(
    file.path(LASAGNA_PROMPTS_DIR, "lasagna_report_data.md")
  )
  text <- omicsai::omicsai_substitute_template(tmpl, c(
    overview,
    list(
      contrast_block        = contrast_block,
      layer_participation   = layer_part,
      modules_summary_table = modsum$table,
      lead_module           = modsum$lead_module,
      module_bridges        = bridges,
      module_detail         = module_detail
    )
  ))

  list(text = text, modules = keep, ctx = ctx)
}

#' Build prompt parameters for a single LASAGNA module summary.
#'
#' Summary mode picks the lead module after community detection — the
#' user-side "contrast" selector chooses the contrast; module choice is
#' driven by the same tier ordering Report mode uses.
lasagna_build_summary_params <- function(res, contrast, pgx, ntop = 12L) {
  ctx <- lasagna_ai_extract_context(res, contrast, pgx, ntop = ntop)
  if (is.null(ctx)) return(NULL)
  ctx$pgx <- pgx

  modules <- .lasagna_compute_modules(res$graph, ctx$nodes_all,
                                      ctx$edges_all,
                                      ntop_nodes = ntop, ntop_edges = ntop)
  if (length(modules) == 0) return(NULL)

  lead_name <- names(modules)[1]
  module_detail <- .lasagna_render_module_block(modules[[lead_name]],
                                                lead_name,
                                                all_modules = modules,
                                                annot = pgx$genes,
                                                contrast = contrast,
                                                ntop = ntop)

  list(
    experiment    = multiomics_ai_experiment_label(pgx),
    contrast      = contrast %||% "N/A",
    module        = lead_name,
    module_detail = module_detail
  )
}

#' Build deterministic methods section for LASAGNA report.
lasagna_build_methods <- function(pgx, contrast) {
  template <- omicsai::omicsai_load_template(
    file.path(LASAGNA_PROMPTS_DIR, "lasagna_methods.md")
  )

  report_date <- format(Sys.Date(), "%Y-%m-%d")
  params <- list(
    experiment = multiomics_ai_experiment_label(pgx),
    contrast   = contrast %||% "N/A",
    date       = report_date
  )

  omicsai::collapse_lines(
    omicsai::omicsai_substitute_template(template, params),
    sprintf("_This report was generated with OmicsPlayground (BigOmics, %s)._",
            report_date),
    "_Note: AI-generated interpretation may contain inaccuracies and must be independently verified._",
    sep = "\n\n"
  )
}
