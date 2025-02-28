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

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
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
      par(mar=c(5,5,1,1))
      ##playbase::wgcna.plotMDS(res, main="", scale=FALSE)
      playbase::wgcna.plotFeatureUMAP(res, nhub=3, method="clust")
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
