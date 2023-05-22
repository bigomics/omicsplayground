##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Expression plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
pcsf_plot_network_ui <- function(id, caption, info.text, height, width) {
  ns <- shiny::NS(id)

  plot_opts <- tagList(
    withTooltip(
      shiny::radioButtons(ns("layout"), "Layout:",
        choiceNames = c("Kamada-Kawai","hierarchical"),
        choiceValues = c("kk","hierarchical"),        
        inline = TRUE),
      "Select layout algorithm",
      placement = "left",
      options = list(container = "body")
    ),
    withTooltip(
      checkboxInput(ns("physics"), "physics", TRUE),
      "....",
      placement = "left",
      options = list(container = "body")
    )
  )

  PlotModuleUI(
    id = ns("plotmodule"),
    title = "PCSF network analysis",
    label = "a",
    plotlib = "visnetwork",
    caption = caption,
    info.text = info.text,
    height = height,
    width = width,
    options = plot_opts,
    download.fmt = c("png", "pdf"),
  )
}

#' PCSF network function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
pcsf_plot_network_server <- function(id,
                                     pgx,
                                     pcsf_compute,
                                     colorby = reactive("gene.cluster"),
                                     contrast = reactive(NULL),
                                     watermark = FALSE
) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    dbg("[pcsf_plot_network.R] ns =",ns("???"))
    
    observeEvent(input$physics, {
      #Update the network
      do.physics <- input$physics
      dbg("[pcsf_plot_network.R] do.physics =",do.physics)      
      visNetwork::visNetworkProxy("plotmodule-renderfigure") %>%
        visNetwork::visPhysics(enabled = do.physics)
    })
    
    get_network <- reactive({
      res <-  pcsf_compute()
      net <- res$net
      .colorby <- colorby()
      .contrast <- contrast()
      if(.colorby=='gene.cluster') {
        igraph::V(net)$type <- igraph::V(net)$group
      } else {
        fx <- res$meta[,.contrast]
        vv <- igraph::V(net)$name
        igraph::V(net)$type <- c("down","up")[1 + 1*(sign(fx[vv])>0)]
      }
      net
    })

    visnetwork.RENDER <- function() {
      net <- get_network()
      dbg("[pcsf_plot_network.R:visnetwork.RENDER] reacted!")

      do.physics <- input$physics
      if(input$layout=="hierarchical") {
        layout <- "hierarchical"
      } else {
        layout <- paste0("layout_with_",input$layout)
      }
      
      ##E(net)$weight <- 1/(E(net)$weight+1e-10)
      visnet <- visplot.PCSF(
        net, style = 1,
        node_size=30, node_label_cex = 30,
        invert.weight = TRUE, edge_width=3,
        Steiner_node_color = "lightblue", Terminal_node_color = "lightgreen",
        extra_node_colors = list("down"="blue", "up"="red"),
        width = '100%', height=900,
        layout = layout,
        physics = do.physics)

      dbg("[pcsf_plot_network.R:visnetwork.RENDER] done!")
      visnet
    }

    
    
    PlotModuleServer(
      "plotmodule",
      func = visnetwork.RENDER,
      plotlib = "visnetwork",
      pdf.width = 10,
      pdf.height = 10,
      add.watermark = watermark,
      vis.delay = 5 ## important! graph physics needs to settle
    )
    
  }) ## end of moduleServer
}
