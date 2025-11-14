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
      inputId = ns("show_all_modules"),
      label = "Show all modules",
      value = FALSE
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
      inputId = ns("topME"),
      label = "Show top ME",
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

      if(n0 > 2) ii <- hclust(Z)$order
      if(n1 > 2) jj <- hclust(tZ)$order

      ## display top 20 modules (or all if <20 available) by default.
      ntop <- if (input$show_all_modules) 10000 else 20
      zm <- rowMeans(sapply(cons$zlist, function(z) rowMeans(z**2)))
      sel <- head(order(-zm), ntop)
      ii <- ii[which(ii %in% sel)]
      zm <- rowMeans(sapply(cons$zlist, function(z) colMeans(z**2)))
      sel <- head(order(-zm),ntop)
      jj <- jj[which(jj %in% sel)]
      
      par(mfrow=c(2,2), mar=c(7,12,2.5,1))
      if(input$transpose) par(mfrow=c(2,2), mar=c(10,9,2.5,1))

      zmax <- max(sapply(cons$zlist, function(z) max(abs(z),na.rm=TRUE)))
      zlim <- c(-zmax, zmax)
      
      ## Individual Module-Trait
      for(i in 1:length(cons$zlist)) {
        k <- names(cons$zlist)[i]
        mat <- cons$zlist[[i]][ii,jj,drop=FALSE]
        if(input$transpose) mat <- Matrix::t(mat)
        main <- paste("module-trait for",toupper(k))
        if(!input$show_all_modules) main <- paste(main, "(top ME)")
        
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
      if(!input$show_all_modules) main <- paste(main, "(top ME)")
      
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
      
      ## ## Discordant Module-Trait
      diffZ <- playbase::wgcna.computeDistinctMatrix(
        cons$zlist,
        ydim = cons$ydim,
        psig = 0.05,
        min.diff = 0.1
      ) 

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
      pdf.width = 10,
      pdf.height = 10,
      res = c(75, 95),
      add.watermark = FALSE
    )


    scatter.RENDER <- function() {
      
      cons <- mwgcna()
      shiny::req(cons)

      F1 <- cor( cons$layers[[1]]$datExpr, cons$layers[[1]]$datTraits, use="pairwise")
      F2 <- cor( cons$layers[[2]]$datExpr, cons$layers[[2]]$datTraits, use="pairwise")
      gg <- names(which(cons$net$colors != "grey"))
      colx <- cons$net$colors
      head(F2)

      kk <- intersect(colnames(F1),colnames(F2))
      if(length(kk)==0) return(NULL)

      kk <- sort(kk)
      F1 <- F1[,kk,drop=FALSE]
      F2 <- F2[,kk,drop=FALSE]      

      F1[is.na(F1)] <- 0
      F2[is.na(F2)] <- 0
      
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
      res = c(90, 110),
      add.watermark = FALSE
    )
    
    
  })
}



