##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_module_graph_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height,
    width) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

wgcna_plot_module_graph_server <- function(id,
                                           wgcna,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    moduleGraph.RENDER <- function() {
      require(igraph)
      res <- wgcna()
      par(mar=c(0,0,0,0))
      playbase::wgcna.plotEigenGeneGraph(res, add_traits=TRUE, main=NULL) 
    }

    PlotModuleServer(
      "plot",
      func = moduleGraph.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 90),
      add.watermark = watermark
    )
  })
}
