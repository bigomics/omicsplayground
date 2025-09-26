##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

consensusWGCNA_plot_moduletrait_ui <- function(
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

consensusWGCNA_plot_moduletrait_server <- function(id,
                                               mwgcna,
                                               r_layers
                                               ) {
  moduleServer(id, function(input, output, session) {
    
    plot.RENDER <- function() {

      cons <- mwgcna()
      shiny::req(cons)
      
      ## fix ordering of heatmaps
      Z <- Reduce('+', lapply(cons$zlist,dist))
      tZ <- Reduce('+', lapply(lapply(cons$zlist,t),dist))
      
      n0 <- nrow(cons$zlist[[1]])
      n1 <- ncol(cons$zlist[[1]])      
      ii <- 1:n0
      jj <- 1:n1      

      if(n0>2) ii <- hclust(Z)$order
      if(n1>2) jj <- hclust(tZ)$order
      
      modtop20 <- input$top20 && length(ii) > 20 
      if(modtop20) {
        zm <- rowMeans(sapply( cons$zlist, function(z) rowMeans(z**2)))
        sel <- head(order(-zm),20)
        ii <- ii[which(ii %in% sel)]
      }
      if(input$top20 && length(jj) > 20) {
        zm <- rowMeans(sapply( cons$zlist, function(z) colMeans(z**2)))
        sel <- head(order(-zm),20)
        jj <- jj[which(jj %in% sel)]
      }

      if(input$transpose) {
        par(mfrow=c(2,3), mar=c(10,9,2.5,1))
      } else {
        par(mfrow=c(2,3), mar=c(7,12,2.5,1))     
      }
      
      ## Individual Module-Trait
      for(i in 1:length(cons$zlist)) {
        k <- names(cons$zlist)[i]
        mat <- cons$zlist[[i]][ii,jj,drop=FALSE]
        if(input$transpose) mat <- Matrix::t(mat)
        main <- paste("module-trait for",toupper(k))
        if(modtop20) main <- paste(main, "(top20)")
        
        playbase::wgcna.plotLabeledCorrelationHeatmap(
          mat,
          nSamples = cons$ydim[i],
          cex.lab = 1.2,
          text = input$showvalues,
          pstar = input$showsig,
          cluster=FALSE, setpar=FALSE,
          main = main
        )
      }
      
      ## Consensus Module-Trait
      consZ <- playbase::wgcna.computeConsensusMatrix(
        cons$zlist,
        ydim = cons$ydim,
        psig = 0.05
      )      
      nsamples <- nrow(cons$datTraits)
      mat <- consZ[ii,jj,drop=FALSE]
      if(input$transpose) mat <- Matrix::t(mat)      
      tt <- paste(toupper(names(cons$zlist)), collapse="+")
      main <- paste("consensus for",tt)
      if(modtop20) main <- paste(main, "(top20)")
      
      playbase::wgcna.plotLabeledCorrelationHeatmap(
        mat,
        nsamples,
        setpar = FALSE,
        text = input$showvalues,
        pstar = input$showsig,
        cluster = FALSE,
        cex.lab = 1.2,
        main = main
      )
      
      ## Discordant Module-Trait
      diffZ <- playbase::wgcna.computeDistinctMatrix(
        cons$zlist,
        ydim = cons$ydim,
        psig = 0.05,
        min.diff = 0.1
      ) 
      for(set in names(diffZ)) {
        mat <- diffZ[[set]][ii,jj,drop=FALSE]
        if(input$transpose) mat <- Matrix::t(mat)
        main <- paste("discordant for",toupper(set))
        if(modtop20) main <- paste(main, "(top20)")
        
        playbase::wgcna.plotLabeledCorrelationHeatmap(
          mat,
          nsamples,
          setpar = FALSE,
          text = FALSE,
          pstar = 1,
          cluster = FALSE,
          cex.lab = 1.2,
          main = main
        )
      }
      
      
    }
    
    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 8,
      pdf.height = 12,
      res = c(90, 105),
      add.watermark = FALSE
    )

    
  })
}



