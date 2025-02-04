##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_gclustering_ui <- function(
    id,
    label,
    title,
    info.text,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  umap.opts <- shiny::tagList(
    shiny::selectInput(ns("clust_method"), "method:", choices = c("tsne2d", "umap2d", "pca2d"))
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    options = umap.opts,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

wgcna_plot_gclustering_server <- function(id,
                                          wgcna.compute,
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    umap.RENDER <- shiny::reactive({
      out <- wgcna.compute()

      method <- "umap2d"
      method <- input$clust_method

      par(mfrow = c(1, 1), mar = c(2, 3, 1, 1))
      me1 <- paste0("ME", out$net$colors)
      pos <- out$clust[[method]]

      playbase::pgx.scatterPlotXY.BASE(pos, var = me1, col = out$me.colors)
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = umap.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 80),
      add.watermark = watermark
    )
  })
}
