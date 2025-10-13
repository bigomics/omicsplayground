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
    ns("heatmap"),
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

consensusWGCNA_plot_moduletrait_scatter_ui <- function(
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
    )
  )

  PlotModuleUI(
    ns("scatter"),
    title = title,
    label = label,
    info.text = info.text,
    # options = options,
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

    ##---------------------------------------------------
    ## heatmap 
    ##---------------------------------------------------    
    
    heatmap.RENDER <- function() {

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

      zmax <- max(sapply(cons$zlist, function(z) max(abs(z),na.rm=TRUE)))
      zlim <- c(-zmax, zmax)
      
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
          main = main,
          zlim = zlim
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
        main = main,
        zlim = zlim        
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
          main = main,
          zlim = zlim          
        )
      }

      ## Effective consensus
      z0 <- consZ
      z0[is.na(z0)] <- 0
      z1 <- diffZ
      for(i in 1:length(z1)) z1[[i]][is.na(z1[[i]])] <- 0
      effZ <- abs(z0) - Reduce('+', lapply(z1,abs))
      playbase::wgcna.plotLabeledCorrelationHeatmap(
        effZ[ii,jj],
        nsamples,
        setpar=FALSE,
        colorpal = playbase::purpleGreyYellow,
        text=FALSE,
        pstar=1,
        cluster=FALSE,
        cex.lab=1.2,
        main = "consensus score"
      )
      
      
    }
    
    PlotModuleServer(
      "heatmap",
      func = heatmap.RENDER,
      pdf.width = 8,
      pdf.height = 12,
      res = c(90, 105),
      add.watermark = FALSE
    )

    ##---------------------------------------------------
    ## scatterplot
    ##---------------------------------------------------    

    scatter.RENDER <- function() {
      
      cons <- mwgcna()
      shiny::req(cons)

      names(cons$layers)
      F1 <- cor( cons$layers[[1]]$datExpr, cons$layers[[1]]$datTraits, use="pairwise")
      F2 <- cor( cons$layers[[2]]$datExpr, cons$layers[[2]]$datTraits, use="pairwise")
      gg <- names(which(cons$net$colors != "grey"))
      colx <- cons$net$colors
      head(F2)
      
      setname <- names(cons$layers)
      par(mfrow=c(2,3), mar=c(5,5,3,1))
      if(ncol(F2)>6) {
        par(mfrow=c(3,4), mar=c(4,5,3,1))
      }
      for(i in 1:ncol(F2)) {
        tt <- colnames(F2)[i] 
        plot( F1[,i], F2[,i], col=colx, main=tt,
          xlab = paste0("trait correlation (",setname[1],")"),
          ylab = paste0("trait correlation (",setname[2],")")
        )
        abline(v=0, h=0, lty=2, lwd=0.6)
      }
      
    }
    
    PlotModuleServer(
      "scatter",
      func = scatter.RENDER,
      pdf.width = 12,
      pdf.height = 8,
      res = c(110, 130),
      add.watermark = FALSE
    )
    
    
  })
}



