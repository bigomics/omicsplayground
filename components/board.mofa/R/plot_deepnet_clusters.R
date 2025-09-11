##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

plot_deepnet_clusters_ui <- function(
    id,
    title = "",
    info.text = "",
    info.methods,
    info.references,
    caption = "",
    label = "",
    height = c("100%", TABLE_HEIGHT_MODAL),
    width = c("auto", "100%")) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    title = title,
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    caption = caption,
    label = label,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

plot_deepnet_clusters_server <- function(id,
                                         net,
                                         update,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot.RENDER <- function(n = 12) {
      update() ## react on updates
      net <- net()
      nsamples <- ncol(net$X[[1]])
      cex <- ifelse(nsamples < 40, 1.4, 1.15)
      cex <- ifelse(nsamples > 100, 0.9, cex)
      par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
      playbase::deep.plotRedux(net,
        method = "pca",
        views = "multi-omics", par = FALSE, cex = cex
      )
    }

    plot.RENDER2 <- function(n = 12) {
      update() ## react on updates
      net <- net()
      nsamples <- ncol(net$X[[1]])
      cex <- ifelse(nsamples < 40, 1.4, 1.15)
      cex <- ifelse(nsamples > 100, 0.9, cex)
      ntypes <- length(net$X)
      nc <- ceiling(sqrt(ntypes))
      nr <- ceiling(ntypes / nc)
      par(mfrow = c(nr, nc), mar = c(4, 4, 2, 1))
      playbase::deep.plotRedux(
        net,
        method = "pca", views = NULL, par = FALSE, cex = cex
      )
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      func2 = plot.RENDER2,
      pdf.width = 10, pdf.height = 10,
      res = c(75, 120),
      add.watermark = watermark
    )
  })
}
