##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_eigengene_heatmap_ui <- function(
    id,
    label,
    title,
    info.text,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(ns("addtraits"),"Add traits", TRUE)
  )
  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    height = height,
    caption = caption,
    width = width,
    download.fmt = c("png", "pdf", "csv")
  )
}

wgcna_plot_eigengene_heatmap_server <- function(id,
                                                wgcna,
                                                watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    csvFunc <- function() {
      res <- wgcna()
      playbase::wgcna.plotEigenGeneAdjacencyHeatmap(
        res, add_traits = input$addtraits, justdata=TRUE)
    }
    
    plot.RENDER <- function() {
      res <- wgcna()
      playbase::wgcna.plotEigenGeneAdjacencyHeatmap(
        res, add_traits = input$addtraits,
        main = "", marx=0.7)
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      csvFunc = csvFunc,
      pdf.width = 8, pdf.height = 6,
      res = c(80, 100),
      add.watermark = watermark
    )
  })
}
