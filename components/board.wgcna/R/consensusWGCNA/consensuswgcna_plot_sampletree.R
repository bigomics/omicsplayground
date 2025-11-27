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

consensusWGCNA_plot_sampletree_server <- function(id,
                                                  mwgcna
                                                  ) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {

      cons <- mwgcna()
      shiny::req(cons)
      
      nsets <- length(cons$datExpr)
      layout.matrix <- matrix(1:(2*nsets), nrow = 2, ncol = nsets)
      layout(layout.matrix, heights=c(1,2), widths=rep(1,nsets))
      what <- if (input$showtraits) "both" else "me"
      
      for(i in 1:nsets) {
        playbase::wgcna.plotConsensusSampleDendroAndColors(
          cons,
          i,
          main = toupper(names(cons$datExpr)[i]),
          what = what,
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



