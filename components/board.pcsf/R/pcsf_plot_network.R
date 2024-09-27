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
    withTooltip(shiny::radioButtons(ns("layout"), "Layout algorithm:",
      choiceNames = c("Barnes-Hut", "Hierarchical", "Kamada-Kawai"),
      choiceValues = c("BH", "hierarchical", "KK" ),
      selected = "KK",
      inline = FALSE
      ),
      "Select graph layout algorithm. Barnes-Hut is a physics-based force-directed layout that is interactive. The Kamada-Kawai layout is based on a physical model of springs but is static. The hierachical layout places nodes as a hierarchical tree."),
    hr(),
    withTooltip(
      radioButtons(
        ns("highlightby"),
        "Highlight labels by:",
        choices = c("centrality", "foldchange" = "prize"),
        selected = "centrality",
        inline = FALSE
      ),
      "Highlight labels by scaling label size with selection."
    ),
    withTooltip(
      shiny::radioButtons(ns("layout"), "Layout algorithm:",
        choiceNames = c("Barnes-Hut", "Kamada-Kawai", "hierarchical"),
        choiceValues = c("BH", "KK", "hierarchical"),
        selected = "",
        inline = FALSE
      ),
      "Select graph layout algorithm. Barnes-Hut is a physics-based force-directed layout that is interactive. The Kamada-Kawai layout is based on a physical model of springs but is static. The hierachical layout places nodes as a hierarchical tree."
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
                                     r_layout = reactive("KK"),
                                     ##                                     highlightby = reactive("none"),
                                     watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$physics, {
      # Update the network
      do.physics <- input$physics
      visNetwork::visNetworkProxy(ns("plotmodule-renderfigure")) %>%
        visNetwork::visPhysics(enabled = do.physics)
    })

    visnetwork.RENDER <- function() {
      sel.layout <- input$layout
      req(sel.layout, input$highlightby)
      
      physics = TRUE
      if (sel.layout == "hierarchical") {
        layout <- "hierarchical"
        physics <- FALSE
      } else if (sel.layout == "KK") {
        layout <- "layout_with_kk"
        physics <- FALSE
      } else {
        ## barnesHut
        layout <- "layout_with_kk"
        physics <- TRUE
      }

      ## compute PCSF
      pcsf <- pcsf_compute()
      
      plt <- playbase::plotPCSF(
        pcsf,
        highlightby = input$highlightby,
        layout = layout,
        physics = physics,
        plotlib = "visnet",
        node_cex = 30,
        label_cex = 30,
        nlabel = -1
      )

      plt
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
