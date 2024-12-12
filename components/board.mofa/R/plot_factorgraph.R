##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_factorgraph_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::checkboxInput(ns("mst"),"Min. spanning tree (MST)", TRUE)
  )

  PlotModuleUI(
    id = ns("module"),
    plotlib = "visnetwork",
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

mofa_plot_factorgraph_server <- function(id,
                                         mofa,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {
      res <- mofa()
      par(mar=c(0,0,0,0))
      gr <- res$graphs$factors
      igraph::add_shape('dot');
      igraph::add_shape('star')
      igraph::V(gr)$label <- igraph::V(gr)$name        
      mst <- input$mst
      ## plot( gr )
      vis <- playbase::mofa.plot_module(
        gr, mst=mst, cex=0.3, plotlib="visnet")
      vis <- vis %>%
        visNetwork::visPhysics(
          barnesHut = list(
            gravitationalConstant = -1000,
            centralGravity = 0.3,
            springLength = 15
          )
        )
      vis
    }
   
    PlotModuleServer(
      id = "module",
      plotlib = "visnetwork",
      # renderFunc = renderUI,
      func = plot.RENDER,
      pdf.width = 8, pdf.height = 8,
      res = c(80, 100),
      add.watermark = watermark
    )

    
  })
}



