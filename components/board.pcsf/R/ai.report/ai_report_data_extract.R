##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# PCSF AI Report — Extraction Helpers
# =============================================================================

pcsf_get_contrasts <- function(pgx) {
  if (is.null(pgx)) return(character(0))
  ct <- playbase::pgx.getContrasts(pgx)
  sort(ct[!grepl("^IA:", ct)])
}

pcsf_infer_datatypes <- function(pgx) {
  if (is.null(pgx$datatype) || pgx$datatype != "multi-omics") return(NULL)
  datatypes <- unique(playbase::mofa.get_prefix(rownames(pgx$X)))
  if (all(c("gx", "px") %in% datatypes)) {
    datatypes <- setdiff(datatypes, c("gx"))
  }
  datatypes
}

pcsf_compute_graph <- function(pgx,
                               contrast,
                               ntop = 750L,
                               as_prize = "fc",
                               add_vhce = TRUE,
                               meta_solution = FALSE) {
  if (is.null(contrast) || !nzchar(contrast)) return(NULL)

  comparisons <- pcsf_get_contrasts(pgx)
  if (!meta_solution && !contrast %in% comparisons) return(NULL)

  graph <- tryCatch(
    playbase::pgx.computePCSF(
      pgx,
      contrast = if (meta_solution) NULL else contrast,
      datatypes = pcsf_infer_datatypes(pgx),
      as_prize = as_prize,
      use_rank = FALSE,
      ntop = as.integer(ntop),
      ncomp = 5,
      beta = 1,
      rm.negedge = TRUE,
      highcor = if (isTRUE(add_vhce)) 0.9 else Inf,
      dir = "both",
      ppi = c("STRING", "GRAPHITE")
    ),
    error = function(e) NULL
  )

  if (is.null(graph) || !inherits(graph, "igraph")) return(NULL)
  graph
}

pcsf_build_centrality_table <- function(pgx,
                                        contrast,
                                        pcsf_graph,
                                        n = 100L) {
  if (is.null(pcsf_graph) || !inherits(pcsf_graph, "igraph")) return(NULL)
  if (is.null(contrast) || !nzchar(contrast)) return(NULL)

  df <- tryCatch(
    playbase::pgx.getPCSFcentrality(
      pgx,
      contrast = contrast,
      pcsf = pcsf_graph,
      level = "gene",
      n = as.integer(n),
      plot = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
  if ("centrality" %in% colnames(df)) df$centrality[is.nan(df$centrality)] <- 0
  df
}

pcsf_build_pathway_overlap <- function(pgx,
                                       contrast,
                                       network_genes,
                                       min_overlap = 3L,
                                       n_top = 10L) {
  if (is.null(pgx$gsetX) || is.null(pgx$GMT) || length(network_genes) < 5) return(NULL)

  gmt_genes <- rownames(pgx$GMT)
  network_genes <- intersect(network_genes, gmt_genes)
  if (length(network_genes) < 5) return(NULL)

  gs_fc <- tryCatch(playbase::pgx.getMetaMatrix(pgx, level = "geneset")$fc,
    error = function(e) NULL
  )
  if (is.null(gs_fc) || !contrast %in% colnames(gs_fc)) return(NULL)

  gs_vals <- gs_fc[, contrast]
  gs_vals <- gs_vals[!is.na(gs_vals)]
  if (length(gs_vals) == 0) return(NULL)

  gmt_matrix <- tryCatch(
    pgx$GMT[network_genes, , drop = FALSE],
    error = function(e) NULL
  )
  if (is.null(gmt_matrix)) return(NULL)

  # Some datasets store GMT in sparse/vector-like forms that can drop dimensions.
  # Coerce defensively so colSums always receives a 2D structure.
  if (is.null(dim(gmt_matrix))) {
    gmt_matrix <- matrix(gmt_matrix, nrow = length(network_genes), byrow = FALSE)
  }
  gmt_matrix <- as.matrix(gmt_matrix)
  if (length(dim(gmt_matrix)) != 2 || ncol(gmt_matrix) == 0) return(NULL)

  overlap_counts <- tryCatch(colSums(gmt_matrix != 0), error = function(e) NULL)
  if (is.null(overlap_counts) || length(overlap_counts) == 0) return(NULL)
  relevant <- names(overlap_counts[overlap_counts >= as.integer(min_overlap)])
  relevant <- intersect(relevant, names(gs_vals))
  if (length(relevant) == 0) return(NULL)

  gs_sub <- gs_vals[relevant]
  gs_sub <- gs_sub[order(-abs(gs_sub))]
  gs_sub <- head(gs_sub, as.integer(n_top))

  data.frame(
    pathway = sub(".*:", "", names(gs_sub)),
    logFC = as.numeric(gs_sub),
    network_genes = as.integer(overlap_counts[names(gs_sub)]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

extract_pcsf_context_data <- function(pgx,
                                      contrast,
                                      pcsf_graph = NULL,
                                      centrality_table = NULL,
                                      n_top = 15L,
                                      pcsf_params = NULL) {
  if (is.null(contrast) || !nzchar(contrast)) return(NULL)

  params <- modifyList(
    list(ntop = 750L, as_prize = "fc", add_vhce = TRUE, meta_solution = FALSE),
    pcsf_params %||% list()
  )

  graph <- pcsf_graph
  if (is.null(graph) || !inherits(graph, "igraph")) {
    graph <- pcsf_compute_graph(
      pgx = pgx,
      contrast = contrast,
      ntop = params$ntop %||% 750L,
      as_prize = params$as_prize %||% "fc",
      add_vhce = isTRUE(params$add_vhce),
      meta_solution = isTRUE(params$meta_solution)
    )
  }
  if (is.null(graph) || !inherits(graph, "igraph")) return(NULL)

  hubs <- centrality_table
  if (is.null(hubs) || !is.data.frame(hubs) || nrow(hubs) == 0) {
    hubs <- pcsf_build_centrality_table(
      pgx = pgx,
      contrast = contrast,
      pcsf_graph = graph,
      n = max(100L, as.integer(n_top))
    )
  }
  if (is.null(hubs) || !is.data.frame(hubs)) hubs <- data.frame()

  if (nrow(hubs) > 0 && "centrality" %in% colnames(hubs)) {
    hubs <- hubs[order(-hubs$centrality), , drop = FALSE]
    hubs <- head(hubs, as.integer(n_top))
  }

  node_names <- igraph::V(graph)$name
  node_types <- igraph::V(graph)$type
  if (is.null(node_types)) node_types <- rep(NA_character_, length(node_names))

  n_nodes <- igraph::vcount(graph)
  n_edges <- igraph::ecount(graph)
  n_terminal <- sum(node_types == "Terminal", na.rm = TRUE)
  n_steiner <- sum(node_types == "Steiner", na.rm = TRUE)
  n_components <- tryCatch(igraph::components(graph)$no, error = function(e) NA_integer_)
  density <- tryCatch(igraph::edge_density(graph), error = function(e) NA_real_)

  fc <- tryCatch(playbase::pgx.getMetaMatrix(pgx, level = "gene")$fc, error = function(e) NULL)
  up_nodes <- 0L
  down_nodes <- 0L
  if (!is.null(fc) && contrast %in% colnames(fc)) {
    fx <- fc[intersect(node_names, rownames(fc)), contrast]
    fx <- fx[!is.na(fx)]
    up_nodes <- sum(fx > 0, na.rm = TRUE)
    down_nodes <- sum(fx < 0, na.rm = TRUE)
  }

  if (nrow(hubs) > 0) {
    sym <- if ("symbol" %in% colnames(hubs)) hubs$symbol else rownames(hubs)
    ty <- node_types[match(sym, node_names)]
    if (!"node_type" %in% colnames(hubs)) hubs$node_type <- ty
  }

  pathways <- pcsf_build_pathway_overlap(
    pgx = pgx,
    contrast = contrast,
    network_genes = node_names,
    min_overlap = 3L,
    n_top = 10L
  )

  caveats <- character(0)
  if (n_nodes < 20) caveats <- c(caveats, "Small inferred network; biological interpretation may be unstable.")
  if (nrow(hubs) < 5) caveats <- c(caveats, "Few hub genes passed ranking filters.")
  if (is.null(pathways) || nrow(pathways) == 0) caveats <- c(caveats, "No strong pathway overlap found for current network genes.")

  list(
    contrast = contrast,
    experiment = pgx$name %||% pgx$description %||% "omics experiment",
    organism = pgx$organism %||% "human",
    network = list(
      n_nodes = n_nodes,
      n_edges = n_edges,
      n_terminal = n_terminal,
      n_steiner = n_steiner,
      n_components = n_components,
      density = density,
      up_nodes = up_nodes,
      down_nodes = down_nodes
    ),
    hubs = hubs,
    pathways = pathways,
    reliability = list(
      n_hubs_reportable = nrow(hubs),
      frac_steiner = if (n_nodes > 0) n_steiner / n_nodes else 0,
      has_pathway_signal = !is.null(pathways) && nrow(pathways) > 0,
      caveats = caveats
    )
  )
}
