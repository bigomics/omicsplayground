##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

plot_deepnet_gradients_ui <- function(
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

  options <- tagList(
    shiny::checkboxInput(ns("onlypositive"), "positive only", TRUE)
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    caption = caption,
    label = label,
    options = options,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

plot_deepnet_gradients_server <- function(id,
                                          net,
                                          update,
                                          conditions,
                                          datatypes,
                                          phenoFC,
                                          type = c("barplot", "scatter")[1],
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot.RENDER <- function(n = 12) {
      update()

      cond <- conditions()
      dtypes <- datatypes()

      shiny::validate(shiny::need(length(cond) > 0, "Please select a condition"))
      shiny::validate(shiny::need(length(dtypes) > 0, "Please select a datatype"))

      ## check if we are too early after a change
      net <- net()
      shiny::req(dtypes %in% names(net$X))
      shiny::req(cond %in% net$labels[[1]])

      grad <- net$get_gradients()[[1]]
      fc <- NULL

      grad <- grad[dtypes]
      grad <- lapply(grad, function(g) g[, cond, drop = FALSE])
      ngrad <- length(grad) * ncol(grad[[1]])

      nc <- min(ceiling(1.2 * sqrt(ngrad)), ngrad)
      nr <- ceiling(ngrad / nc)

      if (type == "barplot") {
        par(mfrow = c(nr, nc), mar = c(8, 4, 2, 1))
        if (nr > 1) par(mar = c(4, 4, 0.5, 1))
        n <- round(36 / nc)
        cex.names <- 0.7 + 0.12 * nc
        playbase::deep.plotMultiOmicsGradients(
          grad,
          n = n, onlypositive = input$onlypositive,
          par = FALSE, cex.names = cex.names
        )
      }

      if (type == "scatter") {
        fc <- phenoFC()
        par(mfrow = c(nr, nc), mar = c(4, 4, 2, 1))
        playbase::deep.plotGradientVSFoldchange(grad, fc = fc, par = FALSE)
      }
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 12,
      pdf.height = 5,
      res = c(75, 110),
      add.watermark = watermark
    )
  })
}
