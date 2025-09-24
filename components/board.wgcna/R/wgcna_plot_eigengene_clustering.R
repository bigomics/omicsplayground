##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_eigengene_clustering_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height,
    width) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(ns("addtraits"), "Add traits", TRUE)
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
    download.fmt = c("png", "pdf", "svg")
  )
}

wgcna_plot_eigengene_clustering_server <- function(id,
                                                   wgcna.compute,
                                                   watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot.RENDER <- function() {
      res <- wgcna.compute()
      ## Plot the relationships among the eigengenes and the trait
      playbase::wgcna.plotEigenGeneClusterDendrogram(
        res,
        add_traits = input$addtraits, main = ""
      )
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(80, 100),
      add.watermark = watermark
    )
  })
}
