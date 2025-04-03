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
#' @param width
#'
#' @export
TimeSeriesBoard.clustering_plot_ui <- function(
    id,
    label = "label",
    title = "title",
    caption = "Module are identified through k-nearest neighbour (KNN) analysis of expression data. <add how>. Line plots depicting time point-specific average expression profiles of all genes within a defined module.",
    info.text = "Module are identified through k-nearest neighbour (KNN) analysis of expression data. <add how>. The number of modules (KNN) range from 2 to 10 and can be set from the 'Advanced option' menu on the right side. Line plots depicting time point-specific average expression profiles of all genes within a defined module. The gene mapped within each module are reported in the table below.",
    info.methods = "info.methods",
    info.references = list(),
    info.extra_link = "extra.link",
    height = c("calc(100vh - 310px)", TABLE_HEIGHT_MODAL),
    width = c("auto", "100%")
    ) {
  ns <- shiny::NS(id)
    
  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    ## plotlib = "plotly",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height
  )
}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#' @param comp
#' @param pgx
#' @param res
#' @param ii
#' @param watermark
#'
#'
#'
#' @export
TimeSeriesBoard.clustering_server <- function(id,
                                              pgx,
                                              data,
                                              timefactor,
                                              watermark = FALSE) {

  moduleServer(id, function(input, output, session) {

    plot_data <- shiny::reactive({
      res <- data()
      res
    })
    
    render_plot <- function() {
      res <- plot_data()
      shiny::req(res)
      
      time <- res$time
      xx <- res$X
      modules <- factor(res$modules)
      timefactor <- timefactor()
      plottype <- ifelse(timefactor, "parcoord", "continuous")
      
      playbase::plotTimeSeries.modules(
        time, xx, modules, main="", legend=FALSE,
        plottype = plottype
      ) 
    }

    PlotModuleServer(
      "plot",
      func = render_plot,
      ##plotlib = "plotly",
      ##csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(90, 120), ## resolution of plots
      pdf.width = 10,
      pdf.height = 10,
      add.watermark = watermark
    )
    
  }) ## end of moduleServer
}
