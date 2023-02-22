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
compare_plot_cum_fc2_ui <- function(id,
                                   height,
                                   width,
                                   label) {
  ns <- shiny::NS(id)
  info_text <- "Barplot showing the cumulative fold changes on dataset 2"

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
compare_plot_cum_fc2_server <- function(id,
                                       inputData,
                                       dataset2,
                                       input.contrast1,
                                       input.contrast2,
                                       cum_fc,
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
          x = factor(rownames(F2),levels =rownames(F2)),
          y = as.numeric(F2)
        ),
        x = "x",
        y = "y",
        yaxistitle = "Cumulative foldchange",
        xaxistitle = "Genes",
        title = "Dataset 2",
        type = "bar",
        plotRawValues = TRUE
      )
            fig
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
