##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

lasagna_get_contrasts <- function(pgx) {
  ct1 <- colnames(pgx$mofa$contrasts %||% matrix(nrow = 0, ncol = 0))
  ct2 <- tryCatch(colnames(playbase::pgx.getMetaMatrix(pgx)$fc), error = function(e) character(0))
  multiomics_ai_drop_interaction_contrasts(intersect(ct1, ct2))
}

lasagna_ai_extract_context <- function(res, contrast, pgx, ntop = 12L) {
  if (is.null(res) || is.null(res$graph) || !inherits(res$graph, "igraph")) return(NULL)
  if (is.null(contrast) || !nzchar(contrast)) return(NULL)

  graph <- res$graph
  vis <- tryCatch(visNetwork::toVisNetworkData(graph), error = function(e) NULL)
  if (is.null(vis) || is.null(vis$nodes) || is.null(vis$edges)) return(NULL)

  nodes <- vis$nodes
  edges <- vis$edges

  if (!"layer" %in% colnames(nodes)) nodes$layer <- "unknown"
  if (!"fc" %in% colnames(nodes)) nodes$fc <- NA_real_
  if (!"rho" %in% colnames(nodes)) nodes$rho <- NA_real_

  node_layer <- setNames(as.character(nodes$layer), as.character(nodes$id))
  from_layer <- node_layer[as.character(edges$from)]
  to_layer <- node_layer[as.character(edges$to)]
  if (!"connection_type" %in% colnames(edges)) {
    edges$connection_type <- ifelse(from_layer == to_layer, "intra", "inter")
  }
  if (!"rho" %in% colnames(edges)) edges$rho <- NA_real_
  if (!"weight" %in% colnames(edges)) edges$weight <- NA_real_

  inter_edge <- grepl("->", as.character(edges$connection_type), fixed = TRUE)
  if (!any(inter_edge)) inter_edge <- from_layer != to_layer

  layer_counts <- sort(table(as.character(nodes$layer)), decreasing = TRUE)

  node_rank <- multiomics_ai_top_idx_by_abs(ntop, nodes$fc, nodes$rho)
  top_nodes <- nodes[node_rank, , drop = FALSE]
  top_nodes <- top_nodes[, intersect(c("label", "layer", "fc", "rho"), colnames(top_nodes)), drop = FALSE]

  edge_rank <- multiomics_ai_top_idx_by_abs(ntop, edges$rho, edges$weight)
  top_edges <- edges[edge_rank, , drop = FALSE]
  top_edges <- top_edges[, intersect(c("from", "to", "rho", "weight", "connection_type"), colnames(top_edges)), drop = FALSE]

  caveats <- character(0)
  if (igraph::vcount(graph) < 30) {
    caveats <- c(caveats, "Small multipartite graph after pruning; interpretation may be unstable.")
  }
  if (sum(inter_edge, na.rm = TRUE) < 5) {
    caveats <- c(caveats, "Limited cross-layer connectivity signal in the current contrast.")
  }

  list(
    contrast = contrast,
    experiment = multiomics_ai_experiment_label(pgx),
    network = list(
      n_nodes = igraph::vcount(graph),
      n_edges = igraph::ecount(graph),
      n_layers = length(unique(nodes$layer)),
      n_inter_edges = sum(inter_edge, na.rm = TRUE),
      n_intra_edges = sum(!inter_edge, na.rm = TRUE)
    ),
    layer_counts = layer_counts,
    top_nodes = top_nodes,
    top_edges = top_edges,
    caveats = caveats
  )
}
