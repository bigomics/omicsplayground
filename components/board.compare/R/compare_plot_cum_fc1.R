##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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
                                    title,
                                    info.text,
                                    label) {
  ns <- shiny::NS(id)

  PlotModuleUI(ns("plot"),
    title = title,
    plotlib = "plotly",
    label = label,
    info.text = info.text,
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
                                        pgx,
                                        labeltype,
                                        # dataset2,
                                        getMatrices,
                                        watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- reactive({
      res <- getMatrices()
      FC <- cbind(res$F1, res$F2)
      FC
    })

    plot.RENDER <- shiny::reactive({
      # shiny::req(pgx$X)
      # shiny::req(dataset2)
      shiny::req(getMatrices())

      # Get the cumulative fold changes for dataset 1
      res <- getMatrices()
      FC <- cbind(res$F1, res$F2)
      F1 <- res$F1
      ii <- head(order(-rowMeans(FC**2, na.rm = TRUE)), 40)
      ii <- ii[order(rowMeans(FC[ii, ], na.rm = TRUE))]
      F1 <- F1[ii, , drop = FALSE]

      # rename_by
      if(!is.null(rownames(F1))){
        rownames(F1) <- make.names(playbase::probe2symbol(rownames(F1), pgx$genes, labeltype(), fill_na = TRUE), unique = TRUE)
      }

      # Prepare input for the plot
      d <- data.frame(
        x = factor(rownames(F1), levels = rownames(F1)),
        y = F1
      )
      ycols <- colnames(d[, 2:ncol(d)])
      fillcolor <- c(RColorBrewer::brewer.pal(6, "Set2"), RColorBrewer::brewer.pal(9, "Set1"))
      # Call the plot function
      suppressWarnings(
        fig <- playbase::pgx.barplot.PLOTLY(
          data = data.frame(
            x = factor(rownames(F1), levels = rownames(F1)),
            y = F1
          ),
          x = "x",
          y = ycols,
          fillcolor = fillcolor,
          yaxistitle = "log2FC",
          xaxistitle = "",
          title = "",
          grouped = FALSE
        )
      )
      return(fig)
    })

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = plot.RENDER,
      csvFunc = plot_data,
      res = c(80, 98), ## resolution of plots
      pdf.width = 10, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
