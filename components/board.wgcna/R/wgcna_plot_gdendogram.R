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
  height,
  width,
  ...
  ) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    caption = caption,
    info.text = info.text,
    height = height,
    width = width,
    download.fmt = c("png", "pdf"),
    ...
  )
}

wgcna_plot_gdendogram_server <- function(id,
                                         wgcna.compute,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    RENDER <- function() {
      res <- wgcna.compute()
      playbase::wgcna.plotDendroAndColors(res, main="")
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
