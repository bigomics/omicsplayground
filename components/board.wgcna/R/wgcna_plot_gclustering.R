##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_gclustering_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height = 400,
    width = 400,
    ...) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::radioButtons(ns("method"), "Method:", choices=c("umap","mds"),
      selected="umap", inline = TRUE),
    shiny::checkboxInput(ns("showhub"), "Show hubgenes", TRUE)    
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
    download.fmt = c("png", "pdf", "svg"),
    ...
  )

}

wgcna_plot_gclustering_server <- function(id,
                                          wgcna,
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    RENDER <- function() {
      res <- wgcna()
      par(mar=c(3.2,3.2,0.8,0.5), mgp=c(2.2,0.8,0))
      if(input$method == "mds") {
        playbase::wgcna.plotMDS(res, main="", scale=FALSE)
      } else {
        playbase::wgcna.plotFeatureUMAP(
          res,
          nhub = ifelse(input$showhub, 3, 0),
          main = "",
          method = "clust"
          #set.par = FALSE
        )
      }
    }

    PlotModuleServer(
      "plot",
      func = RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(80, 110),
      add.watermark = watermark
    )
  })
}
