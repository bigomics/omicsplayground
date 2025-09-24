##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_membership_v_trait_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height,
    width) {
  ns <- shiny::NS(id)
  
  options <- shiny::tagList(
    shiny::checkboxInput(ns("showallmodules"),"Show all modules", FALSE)
  )
  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    height = height,
    caption = caption,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

wgcna_plot_membership_v_trait_server <- function(id,
                                                 wgcna,
                                                 selected_module,
                                                 selected_trait,
                                                 watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    csvFunc <- function() {
      res <- wgcna()
      module <- selected_module()
      trait <- selected_trait()
      shiny::req(!is.null(input$showallmodules))
      shiny::req(res)
      shiny::req(!is.null(module) &&  module!="")
      shiny::req(!is.null(trait) && trait!="")

      module <- ifelse( input$showallmodules, NULL, module)
      
      df <- playbase::wgcna.getGeneStats(
        res,
        module = module,
        trait=trait,
        plot = FALSE
      ) 
      df
    }
    
    RENDER <- function() {
      res <- wgcna()
      module <- selected_module()
      trait <- selected_trait()
      shiny::req(!is.null(input$showallmodules))
      col <- "black"
      if(input$showallmodules) col <- NULL
      shiny::req(res)
      shiny::req(!is.null(module) &&  module!="")
      shiny::req(!is.null(trait) && trait!="")

      module <- ifelse( input$showallmodules, NULL, module)
      
      par(mar=c(0,0,0,0))
      df <- playbase::wgcna.getGeneStats(
        res,
        module = module,
        trait = trait,
        col = col,
        main = "",
        plot = TRUE
      ) 
    }

    PlotModuleServer(
      "plot",
      func = RENDER,
      csvFunc = csvFunc,
      pdf.width = 5, pdf.height = 5,
      res = c(85, 110),
      add.watermark = watermark
    )
  })
}
