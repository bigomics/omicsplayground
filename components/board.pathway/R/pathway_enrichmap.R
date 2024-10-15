##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


compute_enrichmentmap <- function(pgx, qsig = 0.05, ntop = 120, wt = 1, contrast = NULL, plot = FALSE) {
  meta <- playbase::pgx.getMetaMatrix(pgx, level = "geneset")
  #
  F <- meta$fc
  Q <- meta$qv
  X <- pgx$gsetX
  G <- Matrix::t(pgx$GMT)
  table(rownames(G) == rownames(F))
  G <- G[grep("^PATHWAY|^GOBP|^GOMF|^C5|^GO_BP|^GO_MF", rownames(G)), ]

  ## take most significant genesets
  if (!is.null(contrast)) {
    qmin <- Q[, contrast]
    fmax <- abs(F[, contrast])
  } else {
    qmin <- apply(Q, 1, min, na.rm = TRUE)
    fmax <- apply(abs(F), 1, max, na.rm = TRUE)
  }
  sig1 <- (qmin < qsig)
  table(sig1)

  gg <- rownames(Q)[sig1]
  gg.score <- -log10(qmin[gg]) * fmax[gg]
  gg <- gg[order(-gg.score)]
  gg <- intersect(gg, rownames(G))
  gg <- head(gg, ntop)
  length(gg)

  ## create multi-mode correlation matrix
  rho1 <- rho2 <- rho3 <- 1
  if (ncol(F) > 2) rho1 <- cor(t(F[gg, ]))
  rho2 <- cor(t(as.matrix(G[gg, ])))
  #
  rho3 <- cor(t(X[gg, ]))
  rho1 <- pmax(rho1, 0)
  rho2 <- pmax(rho2, 0)
  rho3 <- pmax(rho3, 0)
  R1 <- (rho1 * rho2 * rho3)**0.33
  diag(R1) <- 0
  ii <- which(rowSums(R1, na.rm = TRUE) > 0)
  R1 <- R1[ii, ii]
  dim(R1)

  qmin <- qmin[rownames(R1)]
  F <- F[rownames(R1), , drop = FALSE]
  Q <- Q[rownames(R1), , drop = FALSE]
  G <- G[rownames(R1), , drop = FALSE]

  ## geneset graph
  require(igraph)
  graph <- igraph::graph_from_adjacency_matrix(
    R1,
    weighted = TRUE,
    diag = FALSE, mode = "undirected"
  )
  #
  V(graph)$size <- as.numeric(Matrix::rowSums(G[rownames(R1), ] != 0))
  #

  pos <- layout_with_fr(graph, weights = E(graph)$weight**wt) ## fast!
  #
  #
  rownames(pos) <- V(graph)$name

  if (plot) {
    plot(graph, vertex.size = 4, layout = pos, vertex.label = NA)
  }

  ## do clustering
  cl <- igraph::cluster_louvain(graph)
  nclust <- length(table(cl$membership))
  nclust

  ## determine term frequency from descriptions
  terms <- tapply(V(graph)$name, cl$membership, function(s) s)
  terms <- lapply(terms, function(x) gsub(".*[:]|Homo.*|Mus.*|[(]GO.*|[+]", "", x))
  terms <- lapply(terms, function(x) unlist(strsplit(tolower(x), split = "[ _-]")))
  terms.freq <- lapply(terms, function(x) rev(sort(table(x))))
  terms.freq2 <- table(unlist(lapply(terms.freq, function(x) unique(names(x)))))
  unique.terms <- names(which(terms.freq2 == 1))
  terms.freq <- lapply(terms.freq, function(f) f[names(f) %in% unique.terms])
  cluster.label <- sapply(terms.freq, function(s) paste(head(names(s), 3), collapse = " "))

  ## compute cluster positions
  cluster.pos <- apply(pos, 2, function(x) tapply(x, cl$membership, mean))
  cluster.sd <- apply(pos, 2, function(x) tapply(x, cl$membership, mad))

  ## marker genes for each cluster
  gset <- V(graph)$name
  gs.genes <- tapply(gset, cl$membership, function(s) names(which(Matrix::colMeans(G[s, ] != 0) > 0.2)))
  unique.genes <- names(which(table(unlist(gs.genes)) == 1))
  marker.genes <- lapply(gs.genes, function(gg) intersect(gg, unique.genes))

  ## create cluster graph
  cluster.graph <- igraph::graph_from_adjacency_matrix(
    diag(nclust),
    diag = FALSE, mode = "undirected"
  )
  V(cluster.graph)$name <- cluster.label
  V(cluster.graph)$size <- as.integer(table(cl$membership))

  gx.meta <- playbase::pgx.getMetaMatrix(pgx, level = "gene")

  res <- list(
    graph = graph,
    pos = pos,
    cluster.graph = cluster.graph,
    cluster.pos = cluster.pos,
    cluster.sd = cluster.sd,
    F = F,
    Q = Q,
    gx.meta = gx.meta,
    marker.genes = marker.genes,
    qmin = qmin
  )

  res
}


plot_enrichmentmap <- function(res, contrast = NULL, qsig = 0.05,
                               clustlabel.y = 1.15, nodelabel.y = 1, colorscale = "RdBu",
                               cex = 1, rx = 3, lwd = 1, title.x = 0.5, title.y = 1.08,
                               title = "Enrichment Map", label = TRUE, max.edges = 5,
                               paper_bgcolor = "white", plot_bgcolor = "white") {
  ## transfer parameters
  graph <- res$graph
  pos <- res$pos
  cluster.graph <- res$cluster.graph
  cluster.pos <- res$cluster.pos
  cluster.sd <- res$cluster.sd
  F <- res$F
  Q <- res$Q
  gx.meta <- res$gx.meta
  marker.genes <- res$marker.genes
  numclust <- length(V(cluster.graph))

  my.colors <- c(
    purple = "#cdceeb",
    purple2 = "#afaec9",
    red = "#ed1f24",
    green = "#6abd45",
    green2 = "#309b47",
    blue = "#3d68b1"
  )

  vertex.color <- my.colors["red"]
  cluster.color <- "rgb(255,255,255)"
  cluster.color <- paste0(rainbow(numclust), "11")

  ## create clusters -------------------------------
  cluster_shapes <- vector("list", numclust)
  k <- 1
  for (k in 1:numclust) {
    sx <- rx * cluster.sd[k, 1]
    sy <- rx * cluster.sd[k, 2]
    x0 <- cluster.pos[k, 1] - sx
    x1 <- cluster.pos[k, 1] + sx
    y0 <- cluster.pos[k, 2] - sy
    y1 <- cluster.pos[k, 2] + sy

    cluster_shapes[[k]] <- list(
      type = "circle",
      xref = "x", x0 = x0, x1 = x1,
      yref = "y", y0 = y0, y1 = y1,
      fillcolor = cluster.color[k],
      line = list(
        color = "rgb(170,170,170)",
        width = 4 * lwd
      ),
      #
      layer = "below"
    )
  }

  ## create edges -------------------------------
  edge_shapes <- list()
  if (max.edges > 0) {
    wt <- E(graph)$weight / max(E(graph)$weight)
    es <- data.frame(get.edgelist(graph), weight = wt)
    es <- es[order(-es$weight), ]
    ii <- unlist(tapply(1:nrow(es), es[, 1], function(i) head(i, max.edges)))
    jj <- unlist(tapply(1:nrow(es), es[, 2], function(i) head(i, max.edges)))
    es <- es[union(ii, jj), ]

    for (i in 1:nrow(es)) {
      v0 <- es[i, 1]
      v1 <- es[i, 2]
      wt <- es[i, 3]
      edge <- list(
        type = "line",
        line = list(
          color = my.colors["green"],
          width = 1 + 3 * lwd * wt
        ),
        opacity = 0.85,
        x0 = pos[v0, 1],
        y0 = pos[v0, 2],
        x1 = pos[v1, 1],
        y1 = pos[v1, 2],
        layer = "below"
      )
      edge_shapes[[i]] <- edge
    }
  }

  ## create dataframes for plotting ------------------------
  pt.size <- V(graph)$size
  df <- data.frame(
    x = pos[, 1],
    y = pos[, 2],
    text = V(graph)$name,
    text2 = V(graph)$name,
    size = 4 + 15 * (pt.size / mean(pt.size))**0.5,
    color = my.colors["red"]
  )
  df$size <- cex * df$size

  cl.size <- V(cluster.graph)$size
  cl.size <- 10 + 10 * cl.size**0.5
  cluster.df <- data.frame(
    x = cluster.pos[, 1],
    y = cluster.pos[, 2],
    dy = rx * cluster.sd[, 2],
    ylabel <- cluster.pos[, 2] + clustlabel.y * rx * cluster.sd[, 2],
    text = V(cluster.graph)$name,
    size = 8 * cl.size,
    color = cluster.color
  )

  ## cluster nodes with foldchange
  if (!is.null(contrast)) {
    genesets <- V(graph)$name
    fc <- F[genesets, contrast]
    qv <- Q[genesets, contrast]
    df$color <- -log10(1e-4 + qv) * sign(fc) * abs(fc)**0.5
    df$text2 <- paste0(
      "<b>", df$text, "</b>",
      "<br>Enrichment score: ", round(fc, digits = 4),
      "<br>q-value: ", format(qv, digits = 4)
    )
  } else {
    genesets <- V(graph)$name
    F1 <- F[genesets, ]
    Q1 <- Q[genesets, ]
    ct <- paste(colnames(F1), collapse = ", ")
    ff <- apply(round(F1, 3), 1, paste, collapse = ", ")
    qq <- apply(round(Q1, 3), 1, paste, collapse = ", ")
    df$text2 <- paste0(
      "<b>", df$text, "</b>",
      "<br>Contrasts: ", ct,
      "<br>Enrichment score: ", ff,
      "<br>q-value: ", qq
    )
  }

  ## create plotly graph -------------------------------
  fig <- plotly::plot_ly(
    df,
    type = "scatter",
    mode = "markers",
    text = ~text
  )

  no.axis <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )

  fig <- fig %>%
    plotly::layout(
      #      title = list(



      #      ),
      paper_bgcolor = paper_bgcolor,
      plot_bgcolor = plot_bgcolor,
      margin = list(l = 20, r = 20, b = 20, t = 70),
      xaxis = no.axis,
      yaxis = no.axis,
      shapes = c(
        edge_shapes,
        cluster_shapes
      )
    )

  ## plot title
  midx <- (1 - title.x) * min(pos[, 1]) + (title.x) * max(pos[, 1])
  xalign <- cut(title.x, breaks = c(-99, 0.25, 0.75, 99), labels = c("left", "center", "right"))
  fig <- fig %>% plotly::layout(
    annotations = list(
      x = midx,
      y = title.y,
      text = sprintf("<b>%s</b>", title),
      align = xalign,
      font = list(size = 32 * cex, color = "black"),
      xanchor = "middle",
      yanchor = "top",
      showarrow = FALSE,
      xref = "x",
      yref = "paper"
    )
  )

  ## get top DE genes
  if (!is.null(contrast)) {
    fc.qv <- abs(gx.meta$fc[, contrast] * (gx.meta$qv[, contrast] < qsig))
    fc.qv <- fc.qv[order(-fc.qv)]
    top.genes <- names(which(fc.qv > 0.2))
  } else {
    fc.qv <- abs(gx.meta$fc * (gx.meta$qv < qsig))
    fc.qv <- fc.qv[order(-rowMeans(fc.qv, na.rm = TRUE)), ]
    top.genes <- names(which(apply(fc.qv, 1, max, na.rm = TRUE) > 0.2))
  }
  ## intersection of marker genes and top DE
  marker.genes2 <- lapply(marker.genes, function(g) intersect(top.genes, g))

  ## Add cluster hoverinfo (marker genes)
  marker.genes2 <- lapply(marker.genes2, head, 40)
  marker.genes2 <- sapply(marker.genes2, paste, collapse = " ")
  marker.genes2 <- sapply(marker.genes2, function(s) paste(strwrap(s, 30), collapse = "\n"))
  marker.genes2 <- trimws(marker.genes2)
  i <- 1
  for (i in 1:length(cluster_shapes)) {
    cs <- cluster_shapes[[i]]
    fill0 <- cluster.df$color[i]
    fig <- fig %>%
      plotly::add_polygons(
        x = c(cs$x0, cs$x1, cs$x1, cs$x0),
        y = c(cs$y0, cs$y0, cs$y1, cs$y1),
        line = list(width = 0),
        fillcolor = fill0,
        opacity = 0,
        inherit = FALSE,
        name = paste("<b>MARKER GENES:</b><br>", marker.genes2[i]),
        showlegend = FALSE,
        layer = "below"
      )
  }

  ## draw points (gene sets)
  maxabsf <- max(abs(F))
  fig <- fig %>%
    plotly::add_trace(
      data = df,
      x = ~x,
      y = ~y,
      text = ~text2,
      marker = list(
        color = ~color,
        size = ~size,
        opacity = 0.9,
        line = list(
          color = "rgb(100,100,100)",
          width = 1.3
        ),
        colorscale = colorscale,
        cauto = FALSE,
        cmin = -maxabsf,
        cmax = +maxabsf
      ),
      hovertemplate = "%{text}<extra></extra>",
      showlegend = FALSE
    )

  if (label) {
    df$text3 <- gsub(".*:|_Homo.*|_Mus.*|[(]GO.*", "", df$text)
    fig <- fig %>%
      plotly::add_annotations(
        data = df,
        x = ~x,
        y = ~y,
        yshift = nodelabel.y * 15,
        text = ~text3,
        showarrow = FALSE,
        font = list(size = 12 * cex, color = "black"),
        showlegend = FALSE
      )
  }

  ## Cluster titles
  fig <- fig %>%
    plotly::add_annotations(
      data = cluster.df,
      x = ~x,
      y = ~ylabel,
      text = ~text,
      showarrow = FALSE,
      font = list(size = 28 * cex, color = "black"),
      showlegend = FALSE
    )

  fig
}
