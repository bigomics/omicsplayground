##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

preservationWGCNA_plot_dendrograms_ui <- function(
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

preservationWGCNA_plot_dendrograms_server <- function(id,
                                                      rwgcna) {
  moduleServer(id, function(input, output, session) {
    plot.RENDER <- function() {
      res <- rwgcna()
      shiny::validate(shiny::need(!is.null(res), "Please compute"))

      playbase::wgcna.plotDendroAndTraitCorrelation_multi(
        res$layers,
        show.traits = input$showtraits,
        show.contrasts = input$showcontrasts,
        marAll = c(1, 10, 3, 0.2),
        colorHeightMax = 0.75
      )
      
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 12,
      pdf.height = 8,
      res = c(90, 100),
      add.watermark = FALSE
    )
  })
}
