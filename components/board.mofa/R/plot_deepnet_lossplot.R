##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

plot_deepnet_lossplot_ui <- function(
    id,
    title = "",
    info.text = "",
    info.methods,
    info.references,
    caption = "",
    label = "",
    height = c("100%", TABLE_HEIGHT_MODAL),
    width = c("auto", "100%")) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    title = title,
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    caption = caption,
    label = label,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

plot_deepnet_lossplot_server <- function(id,
                                         net,
                                         update,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot.RENDER <- function() {
      update() ## react on updates
      net <- net()
      n <- length(net$loss_history)
      shiny::validate(shiny::need(n > 0, "please run network"))
      par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
      net$plot_loss()
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 10, pdf.height = 10,
      res = c(75, 100),
      add.watermark = watermark
    )
  })
}
