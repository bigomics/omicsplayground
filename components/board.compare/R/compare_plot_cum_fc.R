##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics Sagl. All rights reserved.
##

#' Expression plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
compare_plot_cum_fc_ui <- function(id,
                                   height,
                                   width) {
  ns <- shiny::NS(id)
  info_text <- ""

  PlotModuleUI(ns("plot"),
    title = "Cumulative foldchange",
    plotlib = "base",
    label = "b",
    info.text = info_text,
    download.fmt = c("png", "pdf", "csv"),
    height = height,
    width = width
  )
}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
compare_plot_cum_fc_server <- function(id,
                                       inputData,
                                       dataset2,
                                       input.contrast1,
                                       input.contrast2,
                                       watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- shiny::reactive({
      ngs1 <- inputData()
      ngs2 <- dataset2()

      ct1 <- head(names(ngs1$gx.meta$meta), 2)
      ct2 <- head(names(ngs2$gx.meta$meta), 2)
      ct1 <- input.contrast1()
      ct2 <- input.contrast2()
      shiny::req(ct1)
      shiny::req(ct2)
      if (!all(ct1 %in% names(ngs1$gx.meta$meta))) {
        return(NULL)
      }
      if (!all(ct2 %in% names(ngs2$gx.meta$meta))) {
        return(NULL)
      }

      F1 <- pgx.getMetaMatrix(ngs1)$fc[, ct1, drop = FALSE]
      F2 <- pgx.getMetaMatrix(ngs2)$fc[, ct2, drop = FALSE]

      gg <- intersect(toupper(rownames(F1)), toupper(rownames(F2)))
      g1 <- rownames(F1)[match(gg, toupper(rownames(F1)))]
      g2 <- rownames(F2)[match(gg, toupper(rownames(F2)))]
      F1 <- F1[g1, , drop = FALSE]
      F2 <- F2[g2, , drop = FALSE]
      colnames(F1) <- paste0("1:", colnames(F1))
      colnames(F2) <- paste0("2:", colnames(F2))

      return(cbind(F1, F2))
    })

    cumfcplot.RENDER <- shiny::reactive({
      F <- plot_data()
      indexes <- substr(colnames(F), 1, 1)
      F1 <- F[, indexes == 1, drop = FALSE]
      F2 <- F[, indexes == 2, drop = FALSE]

      ii <- head(order(-rowMeans(F**2)), 50)
      ii <- ii[order(rowMeans(F[ii, ]))]
      F <- F[ii, , drop = FALSE]
      F1 <- F1[ii, , drop = FALSE]
      F2 <- F2[ii, , drop = FALSE]

      par(mfrow = c(1, 1), mar = c(4.5, 0, 1, 2), mgp = c(2.2, 0.8, 0))
      graphics::layout(matrix(c(1, 2, 3), nrow = 1, byrow = T), widths = c(0.5, 1, 1))

      frame()
      mtext(rownames(F),
        cex = 0.80, side = 2, at = (1:nrow(F) - 0.5) / nrow(F),
        las = 1, line = -12
      )
      col1 <- grey.colors(ncol(F1))
      if (ncol(F1) == 1) col1 <- "grey50"
      pgx.stackedBarplot(F1,
        hz = TRUE, las = 1, col = col1,
        cex.names = 0.01, cex.lab = 1.4, space = 0.25,
        xlab = "cumulative foldchange", ylab = ""
      )
      legend("bottomright", colnames(F1),
        fill = grey.colors(ncol(F1)),
        cex = 0.9, y.intersp = 0.9, inset = c(-0.03, 0.02), xpd = TRUE
      )
      title("DATASET1", line = -0.35, cex.main = 1.2)

      col2 <- grey.colors(ncol(F2))
      if (ncol(F2) == 1) col2 <- "grey50"
      pgx.stackedBarplot(F2,
        hz = TRUE, las = 1, col = col2,
        cex.names = 0.01, cex.lab = 1.4, space = 0.25,
        xlab = "cumulative foldchange", ylab = ""
      )
      legend("bottomright", colnames(F2),
        fill = grey.colors(ncol(F2)),
        cex = 0.9, y.intersp = 0.9, inset = c(-0.03, 0.02), xpd = TRUE
      )
      title("DATASET2", line = -0.35, cex.main = 1.2)
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      plotlib = "base",
      func = cumfcplot.RENDER,
      csvFunc = plot_data,
      res = c(80, 98), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
