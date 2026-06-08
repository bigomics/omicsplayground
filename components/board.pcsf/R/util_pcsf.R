##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' Write an igraph network to a CX2 file (Cytoscape Exchange Format v2)
#'
#' @description Serialize an igraph object to CX2 JSON, the format accepted by
#'   both Cytoscape Web (which does not import GraphML) and Cytoscape Desktop.
#'   Node `x`/`y` attributes, if present, are written as layout coordinates;
#'   all other vertex/edge attributes are emitted with declared datatypes.
#'
#' @param g igraph object
#' @param file output file path
#'
#' @export
pcsf.write_cx2 <- function(g, file) {
  vattr <- igraph::vertex_attr(g)
  eattr <- igraph::edge_attr(g)
  n <- igraph::vcount(g)
  m <- igraph::ecount(g)
  ids <- seq_len(n) - 1L ## 0-based node ids

  cx_type <- function(v) {
    if (is.logical(v)) {
      "boolean"
    } else if (is.integer(v)) {
      "integer"
    } else if (is.numeric(v)) {
      "double"
    } else {
      "string"
    }
  }

  ## x/y become layout coordinates; everything else is a node attribute
  has_xy <- all(c("x", "y") %in% names(vattr))
  node_attrs <- setdiff(names(vattr), c("x", "y"))

  ## igraph layouts span only a few units, but Cytoscape reads x/y as pixels,
  ## so without rescaling all nodes collapse into the centre. Scale the bounding
  ## box (aspect-ratio preserved) to ~75px node spacing, and flip Y because
  ## Cytoscape's Y axis points down (keeps the view matching the app).
  if (has_xy) {
    px <- as.numeric(vattr$x)
    py <- as.numeric(vattr$y)
    span <- max(diff(range(px, na.rm = TRUE)), diff(range(py, na.rm = TRUE)))
    if (!is.finite(span) || span == 0) span <- 1
    scale <- 75 * sqrt(max(n, 1)) / span
    px <- (px - min(px, na.rm = TRUE)) * scale
    py <- (py - min(py, na.rm = TRUE)) * scale
    py <- max(py, na.rm = TRUE) - py
  }

  node_decl <- lapply(node_attrs, function(a) list(d = cx_type(vattr[[a]])))
  names(node_decl) <- node_attrs
  edge_decl <- lapply(names(eattr), function(a) list(d = cx_type(eattr[[a]])))
  names(edge_decl) <- names(eattr)

  nodes <- lapply(seq_len(n), function(i) {
    v <- lapply(node_attrs, function(a) vattr[[a]][i])
    names(v) <- node_attrs
    node <- list(id = ids[i], v = v)
    if (has_xy) {
      node$x <- px[i]
      node$y <- py[i]
    }
    node
  })

  el <- igraph::as_edgelist(g, names = FALSE) ## 1-based row indices
  edges <- lapply(seq_len(m), function(j) {
    v <- lapply(names(eattr), function(a) eattr[[a]][j])
    names(v) <- names(eattr)
    list(id = j - 1L, s = ids[el[j, 1]], t = ids[el[j, 2]], v = v)
  })

  cx <- list(
    list(CXVersion = "2.0", hasFragments = FALSE),
    list(metaData = list(
      list(name = "attributeDeclarations", elementCount = 1L),
      list(name = "nodes", elementCount = n),
      list(name = "edges", elementCount = m)
    )),
    list(attributeDeclarations = list(list(nodes = node_decl, edges = edge_decl))),
    list(nodes = nodes),
    list(edges = edges),
    list(status = list(list(success = TRUE, error = "")))
  )

  jsonlite::write_json(cx, path = file, auto_unbox = TRUE, digits = NA)
}


#' @param x
#'
#' @return
#' @export
visplot.PCSF <- function(
  x, style = 0, edge_width = 5, node_size = 40, node_label_cex = 30,
  Steiner_node_color = "lightblue", Terminal_node_color = "lightgreen",
  Terminal_node_legend = "Terminal", Steiner_node_legend = "Steiner",
  layout = "layout_with_fr", physics = TRUE, layoutMatrix = NULL,
  width = 800, height = 800, invert.weight = FALSE,
  extra_node_colors = list(), ...
) {
  subnet <- x
  if (missing(subnet)) {
    stop("Need to specify the subnetwork obtained from the PCSF algorithm.")
  }
  if (class(subnet)[1] != "PCSF" || class(subnet)[2] != "igraph") {
    stop("The subnetwork must be a \"PCSF\" object derived from an \"igraph\" class.")
  }
  if (edge_width < 1) {
    stop("The edge_width must be greater than 1.")
  }
  if (node_size < 10) {
    stop("The node_size must be greater than 10.")
  }
  prize <- abs(igraph::V(subnet)$prize)
  min1 <- 10
  max1 <- node_size
  r1 <- max1 - min1
  min2 <- min(prize)
  max2 <- max(prize)
  r2 <- max2 - min2
  adjusted_prize <- r1 * (prize - min2) / r2 + min1

  weight <- igraph::E(subnet)$weight
  if (invert.weight) weight <- 1 / (weight + 1e-10)
  min1 <- 1
  max1 <- edge_width
  r1 <- max1 - min1
  min2 <- min(weight)
  max2 <- max(weight)
  r2 <- max2 - min2
  adjusted_weight <- r1 * (weight - min2) / r2 + min1

  ## prepare nodes/edge dataframes
  nodes <- data.frame(1:length(igraph::V(subnet)), igraph::V(subnet)$name)
  names(nodes) <- c("id", "name")
  nodes$group <- igraph::V(subnet)$type
  nodes$size <- adjusted_prize
  nodes$title <- nodes$name
  nodes$label <- nodes$name
  nodes$label.cex <- node_label_cex
  nodes$font.size <- node_label_cex
  edges <- data.frame(igraph::ends(subnet, es = igraph::E(subnet)), adjusted_weight)
  names(edges) <- c("from", "to", "value")
  edges$from <- match(edges$from, nodes$name)
  edges$to <- match(edges$to, nodes$name)

  visNet <- visNetwork::visNetwork(
    nodes, edges,
    width = width, height = height, ...
  ) %>%
    visNetwork::visNodes(shadow = list(enabled = TRUE, size = 12)) %>%
    visNetwork::visEdges(scaling = list(min = 6, max = edge_width * 6))

  if (layout != "hierarchical") {
    visNet <- visNet %>%
      visNetwork::visIgraphLayout(
        layout = layout,
        physics = physics,
        layoutMatrix = layoutMatrix
      )
  }
  if (layout == "hierarchical") {
    visNet <- visNet %>% visNetwork::visHierarchicalLayout(direction = "UD")
  }

  visNet <- visNet %>%
    visNetwork::visPhysics(enabled = physics) %>%
    visNetwork::visOptions(highlightNearest = list(enabled = TRUE)) %>%
    visNetwork::visLegend(width = 0.05, useGroups = TRUE, position = "right")

  visNet
}

plot.pcsf.NOTUSED <- function(net, fx0 = NULL, label.cex = 1) {
  ## take largest connected graph
  csize <- clusters(net)$csize
  csize
  net <- igraph::decompose.graph(net)[[which.max(csize)]]

  c1 <- betweenness(net)
  c2 <- closeness(net)
  c3 <- evcent(net)$vector
  tail(sort(c1))
  tail(sort(c2))
  tail(sort(c3))

  if (is.null(fx0)) fx0 <- igraph::V(net)$prize
  fx0 <- tanh(1.3 * fx0)
  vertex.label.cex <- label.cex

  vertex.size <- igraph::V(net)$prize**0.66
  vertex.color <- "lightblue"
  cpal <- colorRampPalette(c("blue2", "grey90", "red3"))(33)
  vertex.color <- cpal[1 + 16 + round(16 * fx0)]
  edge.width <- (1 - igraph::E(net)$weight / max(igraph::E(net)$weight))

  pos <- igraph::layout_with_graphopt(
    net,
    niter = 5000,
    charge = 0.00001,
    mass = 30,
    spring.length = 1,
    spring.constant = 1
  )

  plot(net,
    vertex.size = vertex.size,
    vertex.color = vertex.color,
    vertex.label.cex = vertex.label.cex,
    vertex.label.dist = 0.16,
    vertex.label.degree = pi / 2,
    vertex.label.family = "sans",
    edge.width = 5 * edge.width,
    layout = pos
  )
}
