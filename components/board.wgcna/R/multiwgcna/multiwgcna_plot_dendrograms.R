##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

multiwgcna_plot_dendrograms_ui <- function(
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
      label = "Show trait correlation",
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

multiwgcna_plot_dendrograms_server <- function(id,
                                               mwgcna,
                                               r_layers
                                               ) {
  moduleServer(id, function(input, output, session) {
    
    plot.RENDER <- function() {
      wgcna <- mwgcna()

      layers <- r_layers()
      shiny::req(layers)
      
      layers <- intersect(layers, names(wgcna))
      wgcna <- wgcna[layers]
      
      shiny::req(length(wgcna)>0)

      nw = length(wgcna)
      onerow=TRUE
      if(onerow) {
        nr=1;nc=nw
      } else {
        nc = ceiling(sqrt(nw))
        nr = ceiling(nw / nc)            
      }

      if(input$showtraits) {
        ## Need to set layout manually
        layout.matrix <- matrix(1:(2*nr*nc), nrow = nr*2, ncol = nc)
        graphics::layout(
          mat = layout.matrix,
          heights = rep(c(1,1),nr), 
          widths  = rep(1, nc)
        ) 
        i=1
        for(i in 1:length(wgcna)) {
          power <- wgcna[[i]]$net$power
          playbase::wgcna.plotDendroAndTraitCorrelation(
            wgcna = wgcna[[i]],
            main = paste0("Dendrogram for ",names(wgcna)[i]," (p=",power,")"),
            show.traits = input$showtraits,
            marAll = c(1,7,1,0),
            #use.tree = input$clusterby,
            setLayout = FALSE
          )  
        }

      } else {
        ## Need to set layout manually
        layout.matrix <- matrix(1:(2*nr*nc), nrow = nr*2, ncol = nc)
        graphics::layout(
          mat = layout.matrix,
          heights = rep(c(2.5,1),nr), 
          widths  = rep(1, nc)
        ) 
        
        i=1
        for(i in 1:length(wgcna)) {
          power <- wgcna[[i]]$net$power
          playbase::wgcna.plotDendroAndColors(
            wgcna = wgcna[[i]],
            main = paste0("Dendrogram for ",names(wgcna)[i]," (p=",power,")"),
            marAll = c(3, 5.5, 3, 1),
            setLayout = FALSE
          )
        }
      }
      
    }
    
    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 12,
      pdf.height = 8,
      res = c(100, 130),
      add.watermark = FALSE
    )

    
  })
}



