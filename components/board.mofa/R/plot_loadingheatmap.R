##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_loadingheatmap_ui <- function(
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

mofa_plot_loadingheatmap_server <- function(id,
                                            mofa,
                                            input_factor = reactive(1),
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {
      res <- mofa()
      shiny::req(res)
      ntop <- 40 / length(res$ww)
      k <- input_factor()
      dbg("[mofa_plot_loading_heatmap_server] k = ", k)
      dbg("[mofa_plot_loading_heatmap_server] names(res) = ", names(res))
      factors <- colnames(res$F)
      dbg("[mofa_plot_loading_heatmap_server] factors = ", factors)
      shiny::req(k %in% factors)
      playbase::mofa.plot_loading_heatmap(
        res, k=k, ntop=ntop, type="splitmap", annot = "scores",
        mar = c(5,5,0,3), annot.ht = 3.5, cexRow=0.9 )
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
