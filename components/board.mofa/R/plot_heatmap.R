##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_heatmap_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
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

mofa_plot_heatmap_server <- function(id,
                                     mofa,
                                     input_factor = reactive(1),
                                     watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {
      res <- mofa()
      ntop <- 40 / length(res$ww)
      k <- as.integer(input_factor())
      playbase::mofa.plot_factor_heatmap(
        res, k=k, ntop=ntop, type="splitmap")
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 8, pdf.height = 8,
      res = c(80, 100),
      add.watermark = watermark
    )

    
  })
}
