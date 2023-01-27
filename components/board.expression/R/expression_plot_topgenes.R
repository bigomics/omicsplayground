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
expression_plot_topgenes_ui <- function(id,
                                        label = "",
                                        height,
                                        width){
  ns <- shiny::NS(id)

  info_text <- "The <strong>Top genes</strong> section shows the average expression plots across the samples for the top differentially (both positively and negatively) expressed genes for the selected comparison from the <code>Contrast</code> settings. Under the plot <i>Settings</i>, users can scale the abundance levels (counts) or ungroup the samples in the plot from the <code>log scale</code> and <code>ungroup samples</code> settings, respectively."

  topgenes_opts <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("gx_logscale"), "log scale", TRUE),
      "Logarithmic scale the counts (abundance levels).",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(shiny::checkboxInput(ns("gx_ungroup"), "ungroup samples", FALSE),
      "Ungroup samples in the plot",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(shiny::checkboxInput(ns("gx_showothers"), "show others", FALSE),
      "Show the 'others' class (if any)",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(ns("pltmod"),
    title = "Expression of top differentially expressed genes",
    label = label,
    plotlib = "base",
    info.text = info_text,
    options = topgenes_opts,
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
#' @param inputData
#' @param res
#' @param ii
#' @param watermark
#'
#'
#'
#' @export
expression_plot_topgenes_server <- function(id,
                                            comp,
                                            inputData,
                                            res,
                                            ii,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    # #calculate required inputs for plotting ---------------------------------

    plot_data <- shiny::reactive({
      comp <- comp() #input$gx_contrast
      ngs <- inputData()
      shiny::req(ngs)

      res <- res()
      if (is.null(res) || nrow(res) == 0) {
        return(NULL)
      }

      ## filter on active rows (using search)
      ## ii  <- genetable$rows_all()
      ii <- ii()
      res <- res[ii, , drop = FALSE]
      if (nrow(res) == 0) {
        return(NULL)
      }


      grouped <- 0
      logscale <- 1
      grouped <- !input$gx_ungroup
      logscale <- input$gx_logscale
      showothers <- input$gx_showothers

      mar1 <- 3.5
      ylab <- ifelse(logscale, "log2CPM", "CPM")

      ny <- nrow(ngs$samples) ## ???!!
      show.names <- ifelse(!grouped & ny > 25, FALSE, TRUE)
      ## nx = ifelse(grouped, ngrp, length(y))
      nx <- ifelse(grouped, 3, ny)
      nc <- 4
      nc <- 8
      if (nx <= 3) nc <- 10
      if (nx > 10) nc <- 5
      if (nx > 25) nc <- 4
      srt <- 35
      sumlen.grpnames <- sum(nchar(strsplit(sub(".*:", "", comp), split = "_vs_")[[1]]))
      if (show.names && sumlen.grpnames <= 20) srt <- 0

      return(list(
        res = res,
        ngs = ngs,
        comp = comp,
        grouped = grouped,
        showothers = showothers,
        ylab = ylab,
        srt = srt,
        logscale = logscale,
        show.names = show.names
      ))
    })


    plotly.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)

      nc <- 8
      mar1 = 3.5

      par(mfrow = c(2, nc), mar = c(mar1, 3.5, 1, 1), mgp = c(2, 0.8, 0), oma = c(0.1, 0.6, 0, 0.6))
      i <- 1
      for (i in 1:nrow(pd[["res"]])) {
        ## if(i > length(top.up)) { frame() }
        ## gene = sub(".*:","",top.up[i])
        gene <- rownames(pd[["res"]])[i]
        pgx.plotExpression(
          pd[["ngs"]],
          pd[["gene"]],
          pd[["comp"]],
          pd[["grouped"]],
          max.points = 200,
          logscale = pd[["logscale"]],
          collapse.others = TRUE,
          showothers = pd[["showothers"]],
          ylab=ylab,
          xlab="",
          srt=pd[["srt"]],
          names = show.names,
          main = ""
        )
        title(gene, cex.main = 1, line = -0.6)
      }
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
      # func2 = modal_plotly.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(90,105), ## resolution of plots
      pdf.width=14, pdf.height=3.5,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}