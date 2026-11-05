##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_module_graph_ui <- function(
  id,
  label = "",
  title = "",
  info.text = "",
  caption = "",
  height,
  width
) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::selectInput(ns("startmodule"), "Start module:", choices=NULL)
  )
  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    options = options,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

wgcna_plot_module_graph_server <- function(id,
                                           wgcna,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    observe({
      res <- wgcna()
      me <- c("<none>",sort(names(res$me.genes)))
      shiny::updateSelectInput(session, "startmodule",
        choices = me, sel = me[1]
      )
    })

    plotRENDER <- function() {
      res <- wgcna()
      start <- input$startmodule
      if(input$startmodule == "<none>") start <- NULL
      par(mar=c(4,4,1,1))
      playbase::wgcna.plotModuleTrajectory(
        wgcna = res, start = start, k = NULL, main="")
    }

    
    PlotModuleServer(
      "plot",
      func = plotRENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 90),
      add.watermark = watermark
    )
  })
}
