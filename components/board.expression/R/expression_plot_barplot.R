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
#' @param heightt
#' @param width
#'
#' @export
expression_plot_barplot_ui <- function(
    id,
    title,
    info.text,
    caption,
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
    ),
    withTooltip(shiny::checkboxInput(ns("barplot_showothers"), "show others", FALSE),
      "Show the 'others' class (if any)",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(ns("pltmod"),
    title = title,
    label = label,
    info.text = info.text,
    plotlib = "plotly",
    caption = caption,
    options = plots_barplot_opts,
    download.fmt = c("png", "pdf"),
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
#' @param sel
#' @param res
#' @param watermark
#'
#'
#'
#' @export
expression_plot_barplot_server <- function(id,
                                           comp,
                                           pgx,
                                           sel,
                                           res,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    # #calculate required inputs for plotting ---------------------------------

    plot_data <- shiny::reactive({
      comp <- comp()
      grouped <- input$barplot_grouped
      logscale <- input$barplot_logscale
      showothers <- input$barplot_showothers
      sel <- sel()
      res <- res()

      shiny::validate(shiny::need(!is.null(sel()), tspan("Please select gene in the table.", js = FALSE)))
      psel <- rownames(res)[sel]
      gene <- psel
      srt <- ifelse(grouped, 0, 35)

      return(list(
        pgx = pgx,
        gene = gene,
        comp = comp,
        sel = sel,
        grouped = grouped,
        logscale = logscale,
        showothers = showothers,
        srt = srt
      ))
    })

    plotly.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)

      if (is.null(pd[["sel"]]) || length(pd[["sel"]]) == 0) {
        frame()
        text(0.5, 0.5, tspan("No gene selected", js = FALSE), col = "black")
        return(NULL)
      }

      if (is.null(res) || is.null(sel)) {
        return(NULL)
      }

      fig <- playbase::pgx.plotExpression(
        pgx = pd[["pgx"]],
        probe = pd[["gene"]],
        comp = pd[["comp"]],
        grouped = pd[["grouped"]],
        showothers = pd[["showothers"]],
        max.points = 200, ## slow!!
        names = TRUE,
        logscale = pd[["logscale"]],
        srt = pd[["srt"]],
        xlab = "",
        plotlib = "plotly"
      )
      fig
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      res = c(80, 95), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
