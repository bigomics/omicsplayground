##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

multiwgcna_plot_modulecorr_ui <- function(
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
      inputId = ns("mergemodules"),
      label = "Merge modules",
      value = FALSE
    ),
    shiny::checkboxInput(
      inputId = ns("showtop"),
      label = "Show top modules (max. 20)",
      value = FALSE
    ),
    shiny::checkboxInput(
      inputId = ns("showvalues"),
      label = "Show correlation values",
      value = FALSE
    ),
    shiny::checkboxInput(
      inputId = ns("addtraits"),
      label = "Add traits",
      value = FALSE
    ),
    shiny::checkboxInput(
      inputId = ns("condition"),
      label = "Condition on phenotype",
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

multiwgcna_plot_modulecorr_server <- function(id,
                                              mwgcna,
                                              r_layers,
                                              r_phenotype
                                              ) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {
      wgcna <- mwgcna()

      ## select layers
      layers <- r_layers()
      phenotype <- r_phenotype()      
      layers <- intersect(layers, names(wgcna))
      wgcna <- wgcna[layers]
      shiny::req(length(wgcna)>0)
      
      if(!input$condition) phenotype <- NULL

      if(input$mergemodules) {
        
        playbase::wgcna.plotEigenGeneAdjacencyHeatmap(
          wgcna,
          multi = TRUE,
          phenotype = phenotype,
          add_traits = input$addtraits,
          main = "Eigengene Clustering",
          nmax = ifelse(input$showtop, 20, -1),
          cex.lab = 0.8,
          cex.text = 0.7,
          mask.intra = FALSE,
          plotDendro = TRUE,
          plotHeatmap = TRUE, 
          colorlabel = TRUE,
          text = input$showvalues,
          pstar = !input$showvalues,
          dendro.horiz = TRUE,
          dendro.width = 0.2,
          dendro.labels = FALSE,
          fixclust = FALSE,
          mar1 = c(6.5, 5, 1.2, 0),
          mar2 = c(8, 12, 3, 2)
        )

      } else {
        
        playbase::wgcna.plotMultiEigengeneCorrelation(
          wgcna,
          addtraits = input$addtraits,
          phenotype = phenotype,
          nmax = ifelse(input$showtop, 20, -1),
          main = NULL,
          showvalues = input$showvalues,
          showsig = !input$showvalues,
          fixcluster = TRUE,
          cex.text = 0.7,
          cex.lab = 0.8,
          setpar = TRUE
        ) 
    
      }
      
    }
    
    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 12,
      pdf.height = 8,
      res = c(95, 100),
      add.watermark = FALSE
    )

    
  })
}



