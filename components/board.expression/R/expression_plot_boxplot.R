##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
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
expression_plot_boxplot_ui <- function(id,
                                       label = "",
                                       height,
                                       width) {
  ns <- shiny::NS(id)
  options <- tagList(
    actionButton(ns("button1"), "some action")
  )

  info_text <- "The top N = {12} differentially (both positively and negatively) expressed gene barplot for the selected comparison under the <code>Contrast</code> settings."

  PlotModuleUI(ns("pltmod"),
    title = "Differential expression",
    label = label,
    plotlib = "base",
    info.text = info_text,
    options = NULL,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height
  )
}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#' @param inputData
#' @param sel
#' @param res
#' @param watermark
#'
#'
#'
#' @export
expression_plot_boxplot_server <- function(id,
                                           inputData,
                                           sel,
                                           res,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    # #calculate required inputs for plotting ---------------------------------

    plot_data <- shiny::reactive({
      ngs <- inputData()
      shiny::req(ngs)

      ## get table
      ## sel=1
      sel <- sel()
      if (is.null(sel) || length(sel) == 0) {
        frame()
        text(0.5, 0.5, "No gene selected", col = "black")
        return(NULL)
      }

      res <- res()
      if (is.null(res) || is.null(sel)) {
        return(NULL)
      }

      psel <- rownames(res)[sel]
      gene <- ngs$genes[1, "gene_name"]
      comp <- 1
      grouped <- TRUE
      logscale <- TRUE
      srt <- 45
      gene <- ngs$genes[psel, "gene_name"]
      comp <- input$gx_contrast
      shiny::req(comp)
      grouped <- input$boxplot_grouped
      logscale <- input$boxplot_logscale
      srt <- ifelse(grouped, 0, 35)

      return(
        ngs = ngs,
        gene = gene,
        comp = comp,
        grouped = grouped,
        logscale = logscale,
        srt = srt
      )
    })

    plotly.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)

      par(mfrow = c(1, 1), mar = c(4, 3, 1.5, 1.5), mgp = c(2, 0.8, 0), oma = c(1, 0.5, 0, 0.5))
      pgx.plotExpression(pd[["ngs"]],
        pd[["gene"]],
        comp = pd[["comp"]],
        grouped = pd[["grouped"]],
        max.points = 200, ## slow!!
        names = TRUE,
        logscale = pd[["logscale"]],
        srt = pd[["srt"]]
      )
    }

    modal_plotly.RENDER <- function() {
      fig <- plotly.RENDER() %>%
        plotly::layout(
          font = list(size = 18),
          legend = list(
            font = list(size = 18)
          )
        )
      fig <- plotly::style(fig, marker.size = 20)
      fig
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "base",
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(80, 95), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
