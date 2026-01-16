##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_sampledendrogram_ui <- function(
  id,
  title = "",
  label = "",
  info.text = "",
  caption = "",
  height,
  width
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
    download.fmt = c("png", "pdf", "csv", "svg")
  )
}

wgcna_plot_sampledendrogram_server <- function(id,
                                               wgcna,
                                               what,
                                               watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    csvFunc <- function() {
      res <- wgcna()
      playbase::wgcna.plotSampleDendroAndColors(
        res,
        show.me = (what=="me"),
        show.traits = (what=="traits"),
        show.contrasts = (what=="contrasts"),        
        justdata = TRUE
      )
    }

    RENDER <- function() {
      res <- wgcna()
      playbase::wgcna.plotSampleDendroAndColors(
        res,
        ##what = what,
        show.me = (what=="me"),
        show.traits = (what=="traits"),
        show.contrasts = (what=="contrasts"),
        main = ""
      )
    }

    PlotModuleServer(
      "plot",
      func = RENDER,
      csvFunc = csvFunc,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 90),
      add.watermark = watermark
    )
  })
}
