##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_module_graph_ui <- function(id, height, width) {
  ns <- shiny::NS(id)

  info_text <- "<b>WGCNA module graph.</b>"

  PlotModuleUI(
    ns("plot"),
    title = "Module graph",
    label = "e",
    info.text = info_text,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

wgcna_plot_module_graph_server <- function(id,
                                           wgcna.compute,
                                           labels2rainbow,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    moduleGraph.RENDER <- shiny::reactive({
      require(igraph)
      out <- wgcna.compute()
      net <- out$net
      datExpr <- out$datExpr
      datTraits <- out$datTraits
      moduleColors <- labels2rainbow(out$net)
      ii <- which(!duplicated(out$net$colors))
      me.colors <- moduleColors[ii]
      names(me.colors) <- out$net$colors[ii]
      me.colors

      ## Recalculate MEs with color as labels
      MEs <- out$net$MEs
      dim(MEs)

      clust <- hclust(dist(t(MEs)))
      clust
      phylo <- ape::as.phylo(clust)
      gr <- igraph::as.igraph(phylo, directed = FALSE)

      is.tip <- grepl("^ME", igraph::V(gr)$name)
      module.nr <- as.integer(sub("^ME|Node.*", "", igraph::V(gr)$name))
      module.size <- table(out$net$colors)
      module.size <- module.size / mean(module.size)
      module.size

      igraph::V(gr)$label <- igraph::V(gr)$name
      igraph::V(gr)$label[!is.tip] <- NA
      igraph::V(gr)$color <- me.colors[as.character(module.nr)]
      igraph::V(gr)$size <- 24 * (module.size[as.character(module.nr)])**0.5
      igraph::V(gr)$size[is.na(igraph::V(gr)$size)] <- 0

      par(mfrow = c(1, 1), mar = c(1, 1, 1, 1) * 0)
      igraph::plot.igraph(
        gr,
        layout = igraph::layout.kamada.kawai,
        vertex.label.cex = 0.8,
        edge.width = 3
      )
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = moduleGraph.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 90),
      add.watermark = watermark
    )
  })
}
