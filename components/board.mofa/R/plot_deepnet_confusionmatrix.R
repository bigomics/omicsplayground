##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

plot_deepnet_confusionmatrix_ui <- function(
  id,
  title = "",
  info.text = "",
  info.methods,
  info.references,
  caption = "",
  label = "",
  height = c("100%", TABLE_HEIGHT_MODAL),
  width = c("auto", "100%")
) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::radioButtons(ns("showset"), "Show set:", c("train", "test"), selected = "train")
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    options = options,
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

plot_deepnet_confusionmatrix_server <- function(id,
                                                net,
                                                update,
                                                watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- reactive({
      update() ## react on updates
      shiny::req(input$showset)
      mat <- playbase::deep.getConfusionMatrix(net(), what = input$showset)
      return(list(matrix = mat[[1]]))
    })

    plot.RENDER <- function(n = 12) {
      pd <- plot_data()
      mat <- pd$matrix
      mr <- min(10, 2 + max(nchar(rownames(mat))) / 1.8)
      par(mfrow = c(1, 1), mar = c(mr, 2, 2, mr))
      klrpal <- RColorBrewer::brewer.pal(n = 8, "Blues")
      playbase::gx.imagemap(t(log(1 + mat)), clust = FALSE, cex = 1.5, col = klrpal)
      idx <- which(!is.na(mat), arr.ind = TRUE)
      text(idx[, 1], idx[, 2], labels = mat[idx], col = "orange", cex = 1.5)
      mtext("actual", side = 3, line = 0.5, cex = 1.5)
      mtext("predicted", side = 2, line = 0.5, cex = 1.5)
    }

    plot.RENDER2 <- function(n = 12) {
      mat1 <- playbase::deep.getConfusionMatrix(net(), what = "train")[[1]]
      mat2 <- playbase::deep.getConfusionMatrix(net(), what = "test")[[1]]

      mr <- min(10, 2 + max(nchar(rownames(mat1))) / 1.8)
      par(mfrow = c(1, 2), mar = 1.8 * c(mr, 2, 2, mr))
      klrpal <- RColorBrewer::brewer.pal(n = 8, "Blues")

      playbase::gx.imagemap(t(log(1 + mat1)), clust = FALSE, cex = 1.5, col = klrpal)
      idx <- which(!is.na(mat1), arr.ind = TRUE)
      text(idx[, 1], idx[, 2], labels = mat1[idx], col = "orange", cex = 1.5)
      mtext("actual", side = 3, line = 0.5, cex = 1.5)
      mtext("predicted", side = 2, line = 0.5, cex = 1.5)

      playbase::gx.imagemap(t(log(1 + mat2)), clust = FALSE, cex = 1.5, col = klrpal)
      idx <- which(!is.na(mat2), arr.ind = TRUE)
      text(idx[, 1], idx[, 2], labels = mat2[idx], col = "orange", cex = 1.5)
      mtext("actual", side = 3, line = 0.5, cex = 1.5)
      mtext("predicted", side = 2, line = 0.5, cex = 1.5)
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      func2 = plot.RENDER2,
      pdf.width = 10, pdf.height = 10,
      res = c(75, 120),
      csvFunc = plot_data,
      add.watermark = watermark
    )
  })
}
