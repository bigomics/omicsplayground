##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_snf_heatmap_ui <- function(
    id,
    title = "",
    info.text = "",
    info.references = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::checkboxInput(ns("split"), "Split by datatype", TRUE)
  )
  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    info.references = info.references,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
  
}

mofa_plot_snf_heatmap_server <- function(id,
                                         mofa,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function(legend=FALSE) {
      res <- mofa()
      snf <- res$snf
      validate(need(!is.null(res), "missing MOFA data."))              
      playbase::snf.heatmap(
        snf, res$X, res$samples, nmax=60, legend=legend,
        do.split = input$split)                 
    }

    plot.RENDER2 <- function() {
      plot.RENDER(TRUE)
    }
    
    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      func2 = plot.RENDER2,      
      pdf.width = 5, pdf.height = 5,
      res = c(80, 100),
      add.watermark = watermark
    )

    
  })
}
