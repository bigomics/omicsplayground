##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_dendrogram_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width  = 400) {

  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(
      inputId = ns("show_cormatrix"),
      label = "Show correlation matrix",
      value = FALSE
    )
  )
  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

mofa_plot_dendrogram_server <- function(id,
                                       mofa,
                                       watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    
    plot.RENDER <- function() {
      res <- mofa()

      shiny::req(res)
      
      TOM <- cor(t(res$X))
      geneTree <- hclust(as.dist(1-TOM), method="complete")
      datatype <- sub(":.*","",rownames(res$W))
      datatypeColors <- rainbow(8)[factor(datatype)]      

      if(input$show_cormatrix) {
        par(mar=c(1,1,0,0))
        WGCNA::TOMplot( TOM**1, geneTree, datatypeColors, main="")
      } else {
        par(mar=c(0,0,0,0))        
        WGCNA::plotDendroAndColors(
          geneTree, datatypeColors, dendroLabels=FALSE,
          groupLabels="Datatype", main=NULL,
          rowText = paste("",datatype),
          rowTextAlignment="center")
      }
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(70, 100),
      add.watermark = watermark
    )
  })
}
