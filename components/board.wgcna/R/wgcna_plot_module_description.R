##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_module_description_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height,
    width) {
  ns <- shiny::NS(id)

  opts <- shiny::tagList(
    shiny::checkboxInput(ns("show_cov"), "covariance", FALSE)
  )

  PlotModuleUI(
    ns("text"),
    outputFunc = htmlOutput,
    title = title,
    label = label,
    info.text = info.text,
    #options = opts,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

wgcna_plot_module_description_server <- function(id,
                                                 wgcna,
                                                 selected_module,
                                                 selected_trait,
                                                 watermark = FALSE) {
  moduleServer(id, function(input, output, session) {


    info.RENDER <- function() {
      wgcna <- wgcna()
      mod <- selected_module()
      trait <- selected_trait()      
      res <- "DESCRIPTION"
      div(shiny::HTML(res), class = "ai-info")
    }

    PlotModuleServer(
      "text",
      func = info.RENDER,
      pdf.width = 8, pdf.height = 5,
      res = c(75, 100),
      add.watermark = watermark
    )
  })
}
