##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

consensusWGCNA_plot_sampletree_ui <- function(
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
      inputId = ns("showtraits"),
      label = "Show traits",
      value = TRUE
    ),
    shiny::checkboxInput(
      inputId = ns("showmodules"),
      label = "Show modules",
      value = TRUE
    ),
    shiny::checkboxInput(
      inputId = ns("split"),
      label = "Split dataset",
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
    download.fmt = c("png", "pdf", "svg")
  )
}

consensusWGCNA_plot_sampletree_server <- function(id,
                                                  mwgcna
                                                  ) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {

      cons <- mwgcna()
      shiny::req(cons)

      ##if(input$top20) {}
      
      nsets <- length(cons$datExpr)
      layout.matrix <- matrix( 1:(2*nsets), nrow = 2, ncol = nsets)
      layout(layout.matrix, heights=c(1,2), widths=rep(1,nsets))
      
      for(i in 1:nsets) {
        dt <- toupper(names(cons$datExpr)[i])
        playbase::wgcna.plotConsensusSampleDendroAndColors(
          cons, i,
          main = toupper(dt),          
          what = "both",
          marAll = c(1.2, 10, 2, 0.3),
          clust.expr = TRUE,
          setLayout = FALSE
        ) 
      }
      
    }
    
    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 12,
      pdf.height = 6,
      res = c(95, 110),
      add.watermark = FALSE
    )

    
  })
}



