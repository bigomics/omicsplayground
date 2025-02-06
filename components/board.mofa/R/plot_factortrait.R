##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_factortrait_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(
      inputId = ns("cluster"),
      label = "Cluster heatmap",
      value = TRUE
    )
  )
  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

mofa_plot_factortrait_server <- function(id,
                                         mofa,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    
    plot.RENDER <- function() {
      shiny::req(mofa())

      par(mar=c(6,5,1,1))    
      playbase::mofa.plot_factor_trait_correlation(
        mofa(),
        main = "",
        par = FALSE,
        cluster = input$cluster,
        type = "wgcna",
        cex.lab=0.9,
        cex_text = NULL )       
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 12,
      pdf.height = 8,
      res = c(75, 110),
      add.watermark = watermark
    )

    
  })
}
