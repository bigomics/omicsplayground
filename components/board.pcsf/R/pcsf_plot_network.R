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
                                     pcsf_beta = reactive(1),
                                     colorby = reactive("gene.cluster"),
                                     contrast = reactive(NULL),
                                     show.centrality = reactive(TRUE),
                                     watermark = FALSE
) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    dbg("[pcsf_plot_network.R] ns =",ns("???"))
    
    observeEvent(input$physics, {
      #Update the network
      do.physics <- input$physics
      visNetwork::visNetworkProxy("plotmodule-renderfigure") %>%
        visNetwork::visPhysics(enabled = do.physics)
    })
    
    get_network <- reactive({

      res <-  pcsf_compute()
      shiny::req(res)
            
      ppi <- res$ppi
      terminals <- res$terminals      
      idx <- res$idx
      
      beta <- as.numeric(pcsf_beta()) - 1      
      net <- PCSF::PCSF(ppi, terminals, w=2, b=exp(beta))
      igraph::V(net)$group <- idx[igraph::V(net)$name]

      ## remove small clusters...
      cmp <- igraph::components(net)
      sel.kk <- which(cmp$csize > 0.10 * max(cmp$csize))
      net <- igraph::subgraph(net, cmp$membership %in% sel.kk)
      class(net) <- c("PCSF","igraph")
      net
    })
    
    visnetwork.RENDER <- function() {

      res <- pcsf_compute()      
      net <- get_network()

      dbg("[pcsf_plot_network.R:visnetwork.RENDER] reacted!")
      dbg("[pcsf_plot_network.R:get_network] 2")

      .colorby <- colorby()
      .contrast <- contrast()
      if(.colorby == 'gene.cluster') {
        igraph::V(net)$type <- igraph::V(net)$group
      } else {
        fx <- res$meta[,.contrast]
        vv <- igraph::V(net)$name
        igraph::V(net)$type <- c("down","up")[1 + 1*(sign(fx[vv])>0)]
      }
      
      do.physics <- input$physics
      if(input$layout=="hierarchical") {
        layout <- "hierarchical"
      } else {
        layout <- paste0("layout_with_",input$layout)
      }
      
      label_cex = 30
      if(show.centrality()) {
        ewt <- 1.0 / igraph::E(net)$weight
        bc <- igraph::page_rank(net, weights=ewt)$vector
        ##bc <- igraph::betweenness(net)
        label_cex <- 30 + 80 * (bc / max(bc))**2
      }
      
      ##E(net)$weight <- 1/(E(net)$weight+1e-10)
      visnet <- visplot.PCSF(
        net, style = 1,
        node_size = 30,
        node_label_cex = label_cex,
        invert.weight = TRUE,
        edge_width = 4,
        Steiner_node_color = "lightblue",
        Terminal_node_color = "lightgreen",
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
