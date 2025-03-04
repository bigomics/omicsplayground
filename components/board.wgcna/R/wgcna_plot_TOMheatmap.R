##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_TOMheatmap_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height = 400,
    width = 400,
    ... ) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plotmodule"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "csv", "svg"),
    ...
  )
}

wgcna_plot_TOMheatmap_server <- function(id,
                                         wgcna.compute,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    csvFunc <- function() {
      playbase::wgcna.plotTOM(res, justdata=TRUE)
    }
    
    RENDER <- function() {
      res <- wgcna.compute()
      playbase::wgcna.plotTOM(res)
    }

    PlotModuleServer(
      "plotmodule",
      func = RENDER,
      csvFunc = csvFunc,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 90),
      add.watermark = watermark
    )
  })
}
