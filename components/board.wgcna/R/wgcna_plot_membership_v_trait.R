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

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    height = height,
    caption = caption,
    width = width,
    download.fmt = c("png", "pdf")
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
      df <- playbase::wgcna.getGeneStats(res, module=module, trait=trait, plot=FALSE) 
      df
    }
    
    RENDER <- function() {
      res <- wgcna()
      module <- selected_module()
      trait <- selected_trait()
      par(mar=c(2,2,1,1))
      df <- playbase::wgcna.getGeneStats(res, module=module, trait=trait,
                                         main="", plot=TRUE) 
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
