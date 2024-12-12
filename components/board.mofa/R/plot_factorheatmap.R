##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_factorheatmap_ui <- function(
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

mofa_plot_factorheatmap_server <- function(id,
                                           mofa,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    
    plot.RENDER <- function() {
      shiny::req(mofa())
      par(mar=c(6,5,2,1))    
      playbase::mofa.plot_factor_correlation(
        mofa(), main=NULL, marx=0.7) 
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 8, pdf.height = 8,
      res = c(70, 110),
      add.watermark = watermark
    )

    
  })
}
