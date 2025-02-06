##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_snf_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- tagList(
  )
  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

mofa_plot_snf_server <- function(id,
                                 mofa,
                                 type = c("affinity","tsne","heatmap")[1],
                                 input_pheno = reactive(1),
                                 watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {
      res <- mofa()
      snf <- res$snf
      validate(need(!is.null(res), "missing MOFA data."))              
      k <- input_pheno()
      shiny::req(k)
      cc <- 'black'
      ph <- factor(input_pheno())
      
      if(type == "affinity") {
        par(mfrow=c(2,2), mar=c(6,1,2,8))
        ndim <- ncol(snf$affinityMatrix[[1]])
        if(ndim>20) par(mar=c(3,1,2,4))
        nmat <- length(snf$affinityMatrix)+1
        if(nmat>4) par(mfrow=c(3,3)) 
        playbase::snf.plot_affinity(snf, k=0.5, par=FALSE) 
      }

      if(type == "tsne") {
        cc <- factor(res$samples[,ph])
        par(mfrow=c(2,2))
        par(mar=c(5,5,2,1))
        for(i in 1:length(snf$posx)) {
          plot( snf$posx[[i]], col=cc, pch=19, cex=1,
               xlab = "TSNE-x", ylab = "TSNE-y" )
          title( names(snf$posx)[i], cex.main=1.4 )
        }
      }

      if(type == "heatmap") {
        playbase::snf.heatmap(snf, res$X, res$samples, nmax=60)                 
      }
      
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(80, 100),
      add.watermark = watermark
    )

    
  })
}
