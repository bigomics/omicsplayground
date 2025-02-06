##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_gsetmofa_traitCor_ui <- function(
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

mofa_plot_gsetmofa_traitCor_server <- function(id,
                                               mofa,
                                               watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {
      mofa <- mofa()      
      shiny::req(mofa)
      gset.mofa <- mofa$gset.mofa      
      par(mar=c(6,5,1,1))
      playbase::mofa.plot_factor_trait_correlation(
        gset.mofa, cex_text=0.6, cex.lab=0.9,
        cluster = input$cluster,
        par=FALSE, main = "")
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 8, pdf.height = 12,
      res = c(75, 110),
      add.watermark = watermark
    )
    
  })
}



