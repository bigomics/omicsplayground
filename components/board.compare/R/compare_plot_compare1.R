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
compare_plot_compare1_ui <- function(id,
                                     height,
                                     title,
                                     info.text,
                                     info.methods,
                                     info.references,
                                     info.extra_link,
                                     width) {
  ns <- shiny::NS(id)

  PlotModuleUI(ns("plot"),
    title = title,
    plotlib = "base",
    label = "a",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    download.fmt = c("png", "pdf", "csv", "svg"),
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
compare_plot_compare1_server <- function(id,
                                         pgx,
                                         contrast1,
                                         hilightgenes,
                                         createPlot,
                                         plottype,
                                         dataset2,
                                         getMatrices,
                                         watermark = FALSE,
                                         labeltype) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    call_createPlot <- function(get_data) {
      ct <- contrast1()
      shiny::req(ct)
      pgx1 <- pgx
      pgx2 <- dataset2()
      all.ct <- names(pgx1$gx.meta$meta)
      ## if (length(ct) == 0) ct <- all.ct[1]
      if (!any(ct %in% all.ct)) {
        return(NULL)
      }
      ct <- intersect(ct, all.ct)
      higenes <- hilightgenes()
      cex.lab <- 1.0
      ntop <- 9999
      type <- plottype()
      mat <- getMatrices()
      target_col <- mat$target_col

      if (length(higenes) <= 3) cex.lab <- 1.3
      data <- createPlot(pgx1, pgx1, pgx2, ct, target_col, type, cex.lab, higenes, ntop, get_data, labeltype)
      data
    }

    plot_data <- shiny::reactive({
      data <- call_createPlot(get_data = TRUE)
      return(data)
    })

    plot.RENDER <- shiny::reactive({
      ## shiny::validate(shiny::need(contrast1(), "Please select contrasts and run 'Compute'"))
      call_createPlot(get_data = FALSE)
    })

    PlotModuleServer(
      "plot",
      plotlib = "base",
      func = plot.RENDER,
      csvFunc = plot_data,
      res = c(90, 110), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
