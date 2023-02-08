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
#' @param heightt
#' @param width
#'
#' @export
expression_plot_barplot_ui <- function(id,
                                       label = "",
                                       height,
                                       width) {
  ns <- shiny::NS(id)

  plots_barplot_opts <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("barplot_grouped"), "grouped", TRUE),
      "Group expression values by conditions.",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(shiny::checkboxInput(ns("barplot_logscale"), "log scale", TRUE),
      "Show logarithmic (log2CPM) expression values.",
      placement = "right", options = list(container = "body")
    )
  )

  info_text <- "The top N = {12} differentially (both positively and negatively) expressed gene barplot for the selected comparison under the <code>Contrast</code> settings."

  PlotModuleUI(ns("pltmod"),
    title = "Differential expression",
    label = label,
    info.text = info_text,
    plotlib = "plotly",
    options = plots_barplot_opts,
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
#' @param comp
#' @param ngs
#' @param sel
#' @param res
#' @param watermark
#'
#'
#'
#' @export
expression_plot_barplot_server <- function(id,
                                           comp,
                                           ngs,
                                           sel,
                                           res,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    # #calculate required inputs for plotting ---------------------------------

    plot_data <- shiny::reactive({
      comp <- comp() # input$gx_contrast
      grouped <- input$barplot_grouped
      logscale <- input$barplot_logscale
      ngs <- ngs()
      sel <- sel()
      res <- res()

      psel <- rownames(res)[sel]
      gene <- ngs$genes[1, "gene_name"]

      gene <- ngs$genes[psel, "gene_name"]
      srt <- ifelse(grouped, 0, 35)

      return(list(
        ngs = ngs,
        gene = gene,
        comp = comp,
        sel = sel,
        grouped = grouped,
        logscale = logscale,
        srt = srt
      ))
    })

    plotly.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)

      if (is.null(pd[["sel"]]) || length(pd[["sel"]]) == 0) {
        frame()
        text(0.5, 0.5, "No gene selected", col = "black")
        return(NULL)
      }

      if (is.null(res) || is.null(sel)) {
        return(NULL)
      }

      par(mfrow = c(1, 1), mar = c(4, 3, 1.5, 1.5), mgp = c(2, 0.8, 0), oma = c(1, 0.5, 0, 0.5))


      fig <- pgx.plotExpression(
        pd[["ngs"]],
        pd[["gene"]],
        comp = pd[["comp"]],
        grouped = pd[["grouped"]],
        max.points = 200, ## slow!!
        names = TRUE,
        logscale = pd[["logscale"]],
        srt = pd[["srt"]]
      )
      fig

    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      # func2 = modal_plotly.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(80, 95), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
