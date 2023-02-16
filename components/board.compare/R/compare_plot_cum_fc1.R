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
compare_plot_cum_fc1_ui <- function(id,
                                   height,
                                   width,
                                   label) {
  ns <- shiny::NS(id)
  info_text <- ""

  PlotModuleUI(ns("plot"),
    title = "Cumulative foldchange",
    plotlib = "plotly",
    label = label,
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
#' @export
compare_plot_cum_fc1_server <- function(id,
                                       inputData,
                                       dataset2,
                                       cum_fc,
                                       input.contrast1,
                                       input.contrast2,
                                       watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    cumfcplot.RENDER <- shiny::reactive({
      F <- cum_fc()
      indexes <- substr(colnames(F), 1, 1)
      F1 <- F[, indexes == 1, drop = FALSE]
      F2 <- F[, indexes == 2, drop = FALSE]

      ii <- head(order(-rowMeans(F**2)), 50)
      ii <- ii[order(rowMeans(F[ii, ]))]
      F <- F[ii, , drop = FALSE]
      F1 <- F1[ii, , drop = FALSE]
      F2 <- F2[ii, , drop = FALSE]

      fig <- pgx.barplot.PLOTLY(
        data = data.frame(
          x = factor(rownames(F1),levels =rownames(F1)),
          y = as.numeric(F1)
        ),
        x = "x",
        y = "y",
        yaxistitle = "Cumulative foldchange",
        xaxistitle = "Genes",
        title = "Dataset 1",
        type = "bar",
        plotRawValues = TRUE
      )

      # fig2 <- pgx.barplot.PLOTLY(
      #   data = data.frame(
      #     x = factor(rownames(F2),levels =rownames(F2)),
      #     y = as.numeric(F2)
      #   ),
      #   x = "x",
      #   y = "y",
      #   yaxistitle = "Cumulative foldchange",
      #   xaxistitle = "Genes",
      #   type = "bar",
      #   plotRawValues = TRUE
      # )

      fig

      # par(mfrow = c(1, 1), mar = c(4.5, 0, 1, 2), mgp = c(2.2, 0.8, 0))
      # graphics::layout(matrix(c(1, 2, 3), nrow = 1, byrow = T), widths = c(0.5, 1, 1))
      #
      # frame()
      # mtext(rownames(F),
      #   cex = 0.80, side = 2, at = (1:nrow(F) - 0.5) / nrow(F),
      #   las = 1, line = -12
      # )
      # col1 <- grey.colors(ncol(F1))
      # if (ncol(F1) == 1) col1 <- "grey50"
      # pgx.stackedBarplot(F1,
      #   hz = TRUE, las = 1, col = col1,
      #   cex.names = 0.01, cex.lab = 1.4, space = 0.25,
      #   xlab = "cumulative foldchange", ylab = ""
      # )
      # legend("bottomright", colnames(F1),
      #   fill = grey.colors(ncol(F1)),
      #   cex = 0.9, y.intersp = 0.9, inset = c(-0.03, 0.02), xpd = TRUE
      # )
      # title("DATASET1", line = -0.35, cex.main = 1.2)
      #
      # col2 <- grey.colors(ncol(F2))
      # if (ncol(F2) == 1) col2 <- "grey50"
      # pgx.stackedBarplot(F2,
      #   hz = TRUE, las = 1, col = col2,
      #   cex.names = 0.01, cex.lab = 1.4, space = 0.25,
      #   xlab = "cumulative foldchange", ylab = ""
      # )
      # legend("bottomright", colnames(F2),
      #   fill = grey.colors(ncol(F2)),
      #   cex = 0.9, y.intersp = 0.9, inset = c(-0.03, 0.02), xpd = TRUE
      # )
      # title("DATASET2", line = -0.35, cex.main = 1.2)
      # p <- grDevices::recordPlot()
      # p
    })

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = cumfcplot.RENDER,
      csvFunc = cum_fc,
      res = c(80, 98), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
