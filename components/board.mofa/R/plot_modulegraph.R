##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_modulegraph_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::checkboxInput(ns("mst"),"Min. spanning tree (MST)", TRUE),
    shiny::checkboxInput(ns("rm_singles"),"Remove singletons", TRUE),
    shiny::checkboxInput(ns("top20"),"Top 20", FALSE)
  )
  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    plotlib = "visnetwork",
    download.fmt = c("png", "pdf")
  )
}

mofa_plot_modulegraph_server <- function(id,
                                         mofa,
                                         input_k = reactive(1),
                                         input_pheno = reactive(NULL),
                                         filter_types = reactive(NULL),
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {
      graphs <- mofa()$graphs
      k <- input_k()
      shiny::req(k %in% names(graphs$features))
      shiny::req(filter_types())
      shiny::validate( shiny::need( length(filter_types())>0,
        "Please select at least one datatype"))

      par(mar=c(0,0,0,0))
      gr <- graphs$features[[k]]
      vtype <- sub(":.*","", igraph::V(gr)$name)
      gr <- igraph::subgraph( gr, vids = vtype %in% filter_types())

      if(input$top20) {
        gr <- igraph::subgraph(gr, vids = head(order(-V(gr)$size),20))
      }
      
      vis <- playbase::mofa.plot_module(
        gr, mst = input$mst,
        rm.single = input$rm_singles,
        nlabel = 100,
        cex = 0.33,
        plotlib = "visnet")
      
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
      "plot",
      func = plot.RENDER,
      plotlib = "visnetwork",
      pdf.width = 8, pdf.height = 8,
      res = c(80, 100),
      add.watermark = watermark
    )

    
  })
}



