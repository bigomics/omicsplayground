##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

multiwgcna_plot_dendrograms_ui <- function(
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
      inputId = ns("showtom"),
      label = "Show TOM heatmap",
      value = TRUE
    ),
    shiny::checkboxInput(
      inputId = ns("showtraits"),
      label = "Show traits",
      value = FALSE
    ),
    shiny::checkboxInput(
      inputId = ns("showcontrasts"),
      label = "Show contrasts",
      value = FALSE
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

multiwgcna_plot_dendrograms_server <- function(id,
                                               mwgcna,
                                               r_layers) {
  moduleServer(id, function(input, output, session) {
    plot.RENDER <- function() {
      wgcna <- mwgcna()$layers
      layers <- r_layers()
      shiny::req(layers)

      layers <- intersect(layers, names(wgcna))
      wgcna <- wgcna[layers]

      shiny::req(length(wgcna) > 0)

      playbase::wgcna.plotMultiDendroAndColors(
        wgcna, 
        show.traits = input$showtraits,
        show.contrasts = input$showcontrasts,
        show.tom = input$showtom,
        show.kme = 0, use.tree = 0,
        colorHeight = 0.5,
        main = names(wgcna),
        marAll = c(1,7,1.5,0),
        cex = 0.7 
      )

      
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 12,
      pdf.height = 8,
      res = c(100, 130),
      add.watermark = FALSE
    )
  })
}
