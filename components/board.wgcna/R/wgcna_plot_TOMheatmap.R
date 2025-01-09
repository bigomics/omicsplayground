##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_TOMheatmap_ui <- function(
    id,
    label,
    title,
    info.text,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

wgcna_plot_TOMheatmap_server <- function(id,
                                         wgcna.compute,
                                         power,
                                         tomtype,
                                         networktype,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    TOMplot.RENDER <- function() {
      shiny::req(power())
      shiny::req(networktype())
      shiny::req(tomtype())      

      out <- wgcna.compute()

      playbase::plotTOMfromResults(
        out,
        power = as.numeric(power()),
        tomtype = tomtype(),
        networktype = networktype(),
        nSelect = 1000
      ) 
    }

    PlotModuleServer(
      "plot",
      func = TOMplot.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 90),
      add.watermark = watermark
    )
  })
}
