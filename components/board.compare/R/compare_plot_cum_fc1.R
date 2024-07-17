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
                                    label) {
  ns <- shiny::NS(id)
  info_text <- "Barplot showing the cumulative fold changes on dataset 1"

  PlotModuleUI(ns("plot"),
    title = "Foldchange (Dataset 1)",
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
                                        pgx,
                                        dataset2,
                                        cum_fc,   
                                        ## compute,
                                        watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## cum_fc_triggered <- shiny::reactiveVal(FALSE)
    ## shiny::observeEvent(compute(), {
    ##   cum_fc_triggered(cum_fc())
    ## })

    cumfcplot.RENDER <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(dataset2)
      ##shiny::req(cum_fc_triggered())
      shiny::req(cum_fc())

      # Get the cumulative fold changes for dataset 1
      ##      FC <- cum_fc_triggered()
      FC <- cum_fc()      
      indexes <- substr(colnames(FC), 1, 1)
      F1 <- FC[, indexes == 1, drop = FALSE]
      ii <- head(order(-rowMeans(FC**2)), 40)
      ii <- ii[order(rowMeans(FC[ii, ]))]
      F1 <- F1[ii, , drop = FALSE]

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
      func = cumfcplot.RENDER,
      csvFunc = cum_fc,
      res = c(80, 98), ## resolution of plots
      pdf.width = 10, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
