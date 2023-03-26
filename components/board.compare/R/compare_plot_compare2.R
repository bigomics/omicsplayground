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
compare_plot_compare2_ui <- function(id,
                                     height,
                                     width) {
  ns <- shiny::NS(id)
  info_text <- "<b>Expression scatter plot.</b> Scatter plots of gene expression."

  PlotModuleUI(ns("plot"),
    title = "Dataset 2",
    plotlib = "base",
    label = "b",
    info.text = NULL,
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
compare_plot_compare2_server <- function(id,
                                         pgx,
                                         input.contrast2,
                                         hilightgenes,
                                         createPlot,
                                         plottype,
                                         dataset2,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- shiny::reactive({
      dataset2()
    })

    scatter2.RENDER <- shiny::reactive({
      pgx1 <- pgx
      pgx2 <- plot_data()
      ct2 <- input.contrast2()
      shiny::req(ct2)
      if (!all(ct2 %in% names(pgx2$gx.meta$meta))) {
        return(NULL)
      }
      type <- plottype()
      higenes <- hilightgenes()
      cex.lab <- 1.0
      ntop <- 9999

      if (length(higenes) <= 3) cex.lab <- 1.3
      createPlot(pgx2, pgx1, pgx2, ct2, type, cex.lab, higenes, ntop)
    })

    PlotModuleServer(
      "plot",
      plotlib = "base",
      func = scatter2.RENDER,
      csvFunc = plot_data,
      res = c(90, 110), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
