##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

preservationWGCNA_plot_overlap_ui <- function(
  id,
  title = "",
  info.text = "",
  caption = "",
  label = "",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  ## options <- shiny::tagList(
  ##   shiny::selectInput(
  ##     inputId = ns("layer"),
  ##     label = "Layer",
  ##     choices = NULL
  ##   )
  ## )

  options <- shiny::tagList(
    shiny::checkboxInput(
      inputId = ns("showdendro"),
      label = "Show dendrogram",
      value = TRUE
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}


preservationWGCNA_plot_overlap_server <- function(id,
                                                  rwgcna) {
  moduleServer(id, function(input, output, session) {
    plot.RENDER <- function() {
      res <- rwgcna()
      ## layer <- r_layer()

      par(mar = c(10, 15, 5, 2))
      playbase::wgcna.plotConsensusOverlapHeatmap(
        res$net, res$layer[[2]]$net,
        ## setLabels = names(pres$layers)[1:2],
        plotDendro = input$showdendro,
        lab.line = c(8, 12),
        setpar = FALSE
      )
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 8,
      pdf.height = 12,
      res = c(90, 110),
      add.watermark = FALSE
    )
  })
}
