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
      inputId = ns("showsig"),
      label = "Show significance",
      value = TRUE
    ),
    shiny::checkboxInput(
      inputId = ns("addtraits"),
      label = "Add traits",
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
                                              r_layers
                                              ) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {
      wgcna <- mwgcna()

      ## select layers
      layers <- r_layers()
      sel.layers <- intersect(layers, names(wgcna))
      wgcna <- wgcna[sel.layers]
      shiny::req(length(wgcna)>0)
      
      if(input$mergemodules) {
        
        playbase::wgcna.plotEigenGeneAdjacencyHeatmap(
          wgcna,
          multi = TRUE,
          add_traits = input$addtraits,
          ##main = NULL,
          method = 2,
          nmax = ifelse(input$showtop, 20, -1),
          mask.intra = FALSE,
          plotDendro = TRUE,
          plotHeatmap = TRUE, 
          colorlabel = TRUE,
          text = input$showvalues,
          pstar = input$showsig,
          dendro.horiz = TRUE,
          dendro.width = 0.2,
          dendro.labels = FALSE,
          mar1 = c(6.5, 5, 1.2, 0),
          mar2 = c(8, 12, 3, 2)
        )

      } else {

        ## Show inter-correlation of modules
        me <- lapply(wgcna, function(w) w$net$MEs)
        comb <- combn(length(me),2)
        ncomb <- ncol(comb)
        ncomb
        nsamples <- nrow(wgcna[[1]]$datExpr)
        Y <- wgcna[[1]]$datTraits
        
        nc <- ceiling(sqrt(ncomb))
        nr <- ceiling(ncomb / nc)
        par(mfrow=c(nr,nc), mar=c(8,10,3,1))
        
        for(k in 1:ncol(comb)) {
          i <- comb[1,k]
          j <- comb[2,k]
          M1 <- me[[i]]
          M2 <- me[[j]]
          if(input$addtraits) {
            M1 <- cbind(M1, Y)
            M2 <- cbind(M2, Y)
          }
          R1 <- cor(M1, M2)
          
          if(input$showtop) {
            NMAX = 20
            ii <- head(order(-apply(abs(R1), 1, max)), NMAX)
            jj <- head(order(-apply(abs(R1), 2, max)), NMAX)          
            R1 <- R1[ii,jj]
          }
          
          playbase::wgcna.plotLabeledCorrelationHeatmap(
            R1,
            nsamples,
            text = input$showvalues,
            pstar = input$showsig,
            cluster = TRUE,
            setpar = FALSE,
            main = paste(names(me)[i],"vs.",names(me)[j])
          )

        }
      }
      
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



