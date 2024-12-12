##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_centrality_ui <- function(
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

mofa_plot_centrality_server <- function(id,
                                        mofa,
                                        input_factor = reactive(NULL),
                                        ## input_pheno = reactive(NULL),
                                        show_types = reactive(NULL),
                                        watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {
      res <- mofa()
      shiny::req(res)
      k <- input_factor()
      shiny::req(k %in% colnames(res$F))
      
      show_types <- show_types()
      y <- res$pheno
      
      par(mar=c(4,4,1,0.5))
      playbase::mofa.plot_centrality(
        res, k=k, y=y, show_types=show_types,
        main = "")
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(70, 120),
      add.watermark = watermark
    )
  })
}
