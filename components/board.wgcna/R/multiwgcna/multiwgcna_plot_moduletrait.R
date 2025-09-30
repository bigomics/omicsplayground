##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

multiwgcna_plot_moduletrait_ui <- function(
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

multiwgcna_plot_moduletrait_server <- function(id,
                                               mwgcna,
                                               r_layers
                                               ) {
  moduleServer(id, function(input, output, session) {
    
    plot.RENDER <- function() {
      wgcna <- mwgcna()

      layers <- r_layers()
      sel.layers <- intersect(layers, names(wgcna))
      wgcna <- wgcna[sel.layers]
      shiny::req(length(wgcna)>0)
      
      tmax = nmax = ifelse( input$showtop, 20, 999)
      plottype = input$plottype
      
      merge_modules = input$mergemodules
      if(merge_modules) {
        par(mfrow=c(1,1))
        if(input$transpose) {
          par(mar=c(14,8,2.5,1))
        } else {
          par(mar=c(8,20,2.5,1))
        }

        playbase::wgcna.plotModuleTraitHeatmap(
          wgcna,
          multi = TRUE,
          setpar = FALSE,
          cluster = TRUE,
          main = paste(names(wgcna),collapse=" + "),
          transpose = !input$transpose,
          colorlabel = TRUE,
          nmax = nmax,
          tmax = nmax,          
          text = input$showvalues,
          pstar = !input$showvalues
        ) 

      } else {

        nc <- ceiling(sqrt(length(wgcna)))
        nr <- ceiling(length(wgcna) / nc)
        par(mfrow=c(nr,nc))
        if(input$transpose) {
          par(mar=c(14,8,2.5,1))
        } else {
          par(mar=c(8,14,2.5,1))
        }
        i=1
        for(i in 1:length(wgcna)) {
          
          playbase::wgcna.plotModuleTraitHeatmap(
            wgcna[[i]],
            multi = FALSE,
            setpar = FALSE,
            cluster = TRUE,
            main = names(wgcna)[i],          
            transpose = !input$transpose,
            colorlabel = TRUE,
            nmax = nmax,
            tmax = nmax,            
            text = input$showvalues,
            pstar = !input$showvalues
          ) 
        }

      }
    }
    
    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 8,
      pdf.height = 12,
      res = c(95, 100),
      add.watermark = FALSE
    )

    
  })
}



