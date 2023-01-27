#' ##
#' ## This file is part of the Omics Playground project.
#' ## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
#' ##
#'
#' #' Expression plot UI input function
#' #'
#' #' @description A shiny Module for plotting (UI code).
#' #'
#' #' @param id
#' #' @param label
#' #' @param height
#' #'
#' #' @export
#' expression_plot_FnName_ui <- function(id,
#'                                       label='',
#'                                       height,
#'                                       width) {
#'   ns <- shiny::NS(id)
#'
#'   info_text = ""
#'
#'   PlotModuleUI(ns(""),
#'                title = "",
#'                label = label,
#'                plotlib = "plotly",
#'                info.text = info_text,
#'                options = NULL,
#'                download.fmt=c("png","pdf","csv"),
#'                height = height,
#'                width = width)
#' }
#'
#' #' Expression plot Server function
#' #'
#' #' @description A shiny Module for plotting (server code).
#' #'
#' #' @param id
#' #'
#' #' @return
#' #' @export
#' expression_plot_FnName_server <- function(id, watermark = FALSE)
#' {
#'   moduleServer( id, function(input, output, session) {
#'
#'
#'         #reactive function listening for changes in input
#'         plot_data <- shiny::reactive({
#'           #code here
#'         })
#'
#'         plot.RENDER <- function() {
#'           pd <- plot_data()
#'           shiny::req(pd)
#'
#'           #plot code here
#'         }
#'
#'         plotly.RENDER <- function() {
#'           pd <- plot_data()
#'           shiny::req(pd)
#'
#'           df <- pd
#'
#'           ## plot as regular plot
#'           plotly::plot_ly(data = df,
#'                           type = '',
#'                           x = "",
#'                           y = "",
#'                           ## hoverinfo = "text",
#'                           hovertext = ~annot,
#'                           marker = list(color = ~color)
#'           )
#'         }
#'
#'         modal_plotly.RENDER <- function() {
#'           plotly.RENDER() %>%
#'             plotly::layout(
#'               ## showlegend = TRUE,
#'               font = list(
#'                 size = 16
#'               )
#'             )
#'         }
#'
#'
#'         PlotModuleServer(
#'           "plot",
#'           plotlib = "plotly",
#'           func = plotly.RENDER,
#'           func2 = modal_plotly.RENDER,
#'           csvFunc = plot_data,   ##  *** downloadable data as CSV
#'           res = c(80,170),                ## resolution of plots
#'           pdf.width = 6, pdf.height = 6,
#'           add.watermark = watermark
#'         )
#'     }## end of moduleServer
#' }
