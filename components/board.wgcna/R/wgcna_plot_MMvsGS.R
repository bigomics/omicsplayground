##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_MMvsGS_ui <- function(
  id,
  label = "",
  title = "",
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
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    ## plotlib = "ggiraph",
    download.fmt = c("png", "pdf", "svg")
  )
}

wgcna_plot_MMvsGS_server <- function(id,
                                     wgcna.compute,
                                     watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    RENDER <- function() {
      res <- wgcna.compute()
      trait <- names(res$datTraits)[1]
      module <- names(res$me.colors)[2]
      playbase::wgcna.plotMMvsGS(res, module, trait, abs = FALSE, plotlib = "ggplot")
    }

    PlotModuleServer(
      "plot",
      ## plotlib = "ggiraph",
      func = RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(90, 105),
      add.watermark = watermark
    )
  })
}
