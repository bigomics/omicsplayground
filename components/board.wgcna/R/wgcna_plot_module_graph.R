##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_module_graph_ui <- function(
  id,
  label = "",
  title = "",
  info.text = "",
  caption = "",
  height,
  width
) {
  ns <- shiny::NS(id)

  ## The radio and the conditional panels share this module's namespace;
  ## conditionalPanel() adds that namespace to the input ids in its
  ## condition itself.
  options <- shiny::tagList(
    shiny::radioButtons(ns("plottype"), "Plot type:", c("graph", "TSP"),
      inline = TRUE
    ),
    shiny::conditionalPanel(
      "input.plottype == 'graph'",
      ns = ns,
      shiny::checkboxInput(ns("addtraits"), "Add traits", TRUE)
    ),
    shiny::conditionalPanel(
      "input.plottype == 'TSP'",
      ns = ns,
      shiny::selectInput(
        ns("startmodule"), "Start module:",
        choices = stats::setNames("<none>", "Start from first module")
      )
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    options = options,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

wgcna_plot_module_graph_server <- function(id,
                                           wgcna,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    ## <none> means "start the trajectory at the largest module".
    start_module_choices <- function(me) {
      stats::setNames(me, c("Start from first module", me[-1]))
    }

    ## Start-module choices for the trajectory view. Only invalidate the
    ## select's value when the dataset actually changes it: an unconditional
    ## update would re-render the plot every time wgcna() is re-evaluated.
    shiny::observe({
      res <- wgcna()
      modules <- sort(names(res$me.genes))
      me <- c("<none>", modules)
      current <- input$startmodule
      if (identical(current, me[1]) || (!is.null(current) && current %in% me)) {
        return()
      }
      shiny::updateSelectInput(session, "startmodule",
        choices = start_module_choices(me), selected = me[1]
      )
    })

    plotRENDER <- function() {
      res <- wgcna()
      if (identical(input$plottype, "TSP")) {
        ## TSP
        start <- input$startmodule
        if (is.null(start) || !nzchar(start) || identical(start, "<none>")) start <- NULL
        par(mar = c(4, 4, 1, 1))
        playbase::wgcna.plotModuleTrajectory(
          wgcna = res, start = start, k = NULL, main = ""
        )
      } else {
        par(mar = c(0, 0, 0, 0))
        playbase::wgcna.plotEigenGeneGraph(res, add_traits = input$addtraits,
          main = NULL
        )
      }
    }

    PlotModuleServer(
      "plot",
      func = plotRENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 90),
      add.watermark = watermark
    )
  })
}
