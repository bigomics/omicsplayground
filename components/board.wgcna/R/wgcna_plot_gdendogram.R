##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_gdendogram_ui <- function(
  id,
  title = "",
  label = "",
  info.text = "",
  caption = "",
  height = 400,
  width = 400,
  ...
) {
  ns <- shiny::NS(id)


  options <- shiny::tagList(
    shiny::checkboxInput(ns("showtrait"), "Show traits", FALSE),
    shiny::checkboxInput(ns("showcontrasts"), "Show contrasts", FALSE)
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    options = options,
    label = label,
    caption = caption,
    info.text = info.text,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg"),
    ...
  )
}

wgcna_plot_gdendogram_server <- function(id,
                                         wgcna.compute,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    RENDER <- function() {
      res <- wgcna.compute()
      playbase::wgcna.plotDendroAndColors(
        res,
        show.traits = input$showtrait,
        show.contrasts = input$showcontrasts,
        marAll = c(0.4, 5, 1, 0.2),
        main = "")
    }

    PlotModuleServer(
      "plot",
      func = RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(85, 120),
      add.watermark = watermark
    )
  })
}
