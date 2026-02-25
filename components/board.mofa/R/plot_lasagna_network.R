##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_lasagna_network_ui <- function(
  id,
  title = "",
  info.text = "",
  info.methods = "",
  info.references = NULL,
  info.extra_link = NULL,
  caption = info.text,
  label = "",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  options <- tagList()

  PlotModuleUI(
    ns("plotmodule"),
    plotlib = "visnetwork",
    title = title,
    label = "",
    caption = caption,
    info.text = info.text,
    info.references = info.references,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}


mofa_plot_lasagna_network_server <- function(id,
                                             data,
                                             pgx,
                                             watermark = FALSE) {

  moduleServer(id, function(input, output, session) {
  
    plot.RENDER <- function() {
      res <- data()
      shiny::req(res)
      graph <- res$graph
      dbg("[mofa_plot_lasagna_network_server] layers=", graph$layers)
      min_rho <- 0.0
      ewt <- igraph::E(graph)$weight
      graph <- igraph::subgraph_from_edges(graph, which(abs(ewt) > 0))

      vis <- playbase::lasagna.plot_visgraph(
        graph,
        layers = NULL,
        ntop = -1,
        ecex = 1,
        vcex = 1,
        min_rho = min_rho,
        mst = TRUE
      )

      return(vis)
    }

    PlotModuleServer(
      "plotmodule",
      func = plot.RENDER,
      plotlib = "visnetwork",
      pdf.width = 12, pdf.height = 6,
      res = c(75, 90),
      add.watermark = watermark
    )
  })
}
