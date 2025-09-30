##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

preservationWGCNA_plot_moduletrait_ui <- function(
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
      inputId = ns("top20"),
      label = "Show top 20",
      value = TRUE
    ),
    shiny::checkboxInput(
      inputId = ns("showvalues"),
      label = "Show correlation values",
      value = FALSE
    ),
    shiny::checkboxInput(
      inputId = ns("showsig"),
      label = "Show p-values",
      value = TRUE
    ),
    shiny::checkboxInput(
      inputId = ns("transpose"),
      label = "Transpose matrix",
      value = FALSE
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
    download.fmt = c("png", "pdf", "svg")
  )
}

preservationWGCNA_plot_moduletrait_server <- function(id,
                                               rwgcna
                                               ) {
  moduleServer(id, function(input, output, session) {
    
    plot.RENDER <- function() {

      res <- rwgcna()
      shiny::req(res)

      par(mfrow=c(2,2), mar=c(10,9,3,2))
      playbase::wgcna.plotPreservationModuleTraits(
        res,
        subplots = FALSE,
        setpar = FALSE,
        rm.na = TRUE
      ) 
      
      
    }
    
    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 8,
      pdf.height = 12,
      res = c(80, 100),
      add.watermark = FALSE
    )

    
  })
}



