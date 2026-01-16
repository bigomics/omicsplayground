##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

consensusWGCNA_plot_dendrograms_ui <- function(
  id,
  title = "",
  info.text = "",
  caption = "",
  label = "",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(
      inputId = ns("showtraits"),
      label = "Show traits",
      value = TRUE
    ),
    shiny::checkboxInput(
      inputId = ns("showcontrasts"),
      label = "Show contrasts",
      value = TRUE
    ),
    shiny::selectInput(
      inputId = ns("clusterby"),
      label = "Cluster by",
      choices = NULL
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

consensusWGCNA_plot_dendrograms_server <- function(id,
                                                   mwgcna) {
  moduleServer(id, function(input, output, session) {
    observeEvent(mwgcna(), {
      cons <- mwgcna()
      shiny::req(cons)
      trees <- c("Consensus" = 0, names(cons$layers))
      shiny::updateSelectInput(session, "clusterby", choices = trees)
    })

    plot.RENDER <- function() {
      cons <- mwgcna()
      shiny::req(cons)

      mytrees <- c(0, names(cons$layers))
      shiny::req(input$clusterby %in% mytrees)

      playbase::wgcna.plotDendroAndTraitCorrelation(
        cons,
        main = "",
        show.traits = input$showtraits,
        show.contrasts = input$showcontrasts,
        marAll = c(1, 6, 1, 0),
        use.tree = input$clusterby,
        colorHeightMax = 0.75
      )
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 12,
      pdf.height = 8,
      res = c(90, 110),
      add.watermark = FALSE
    )
  })
}
