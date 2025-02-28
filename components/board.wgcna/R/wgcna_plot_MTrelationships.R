##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_MTrelationships_ui <- function(
    id,
    title = "",
    label = "",
    info.text = "",
    caption = "",
    height,
    width) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(ns("cluster"), "Cluster heatmap", TRUE)
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
    # FIXME png and pdf is not working, to avoid crash, we decided to remove it    
    download.fmt = c("png","pdf","csv","svg")
  )
}

wgcna_plot_MTrelationships_server <- function(id,
                                              wgcna.compute,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    csvFunc <- function() {
      res <- wgcna.compute()
      playbase::wgcna.plotModuleTraitHeatmap(
        res, setpar=FALSE, cluster=input$cluster,
        justdata = TRUE) 
    }
    
    RENDER <- function() {
      res <- wgcna.compute()
      par(mar = c(6, 8, 1, 0.4))
      playbase::wgcna.plotModuleTraitHeatmap(
        res, setpar=FALSE, cluster=input$cluster, main='') 
    }

    PlotModuleServer(
      "plot",
      func = RENDER,
      csvFunc = csvFunc,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 100),
      add.watermark = watermark
    )
  })
}
