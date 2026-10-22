##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
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
