##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_clustering_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

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

mofa_plot_clustering_server <- function(id,
                                        data,
                                        type = c("samples","features")[1],
                                        input_contrast = reactive(NULL),
                                        watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {
      res <- data()
      k <- input_contrast()
      shiny::req(res$posx)
      shiny::req(res$posf)
      shiny::req(res$Y)      

      dbg("[mofa_clustering_server] k = ", k)      
      dbg("[mofa_clustering_server] colnames.Y = ", colnames(res$Y))
      dbg("[mofa_clustering_server] k in colnames.Y = ", k %in% colnames(res$Y))
      if(!is.null(k)) shiny::req(k %in% colnames(res$Y))

      dbg("[mofa_clustering_server] dim(res$posx) = ", dim(res$posx))
      dbg("[mofa_clustering_server] dim(res$posf) = ", dim(res$posf))      
      
      col1 = 'black'      
      if( type == "samples") {
        shiny::req(res$posx)
        par(mfrow=c(2,2), mar=c(4,4,2.5,1))
        if(length(res$posx)>4) par(mfrow=c(3,3))
        for(i in 1:length(res$posx)) {
          plot( res$posx[[i]], col = col1, pch=20, cex=2.2,
               xlab="UMAP1", ylab="UMAP2" )
          title(toupper(names(res$posx)[i]), cex.main=1.2)
        }
      }
      if( type == "features") {

        ## color by contrast correlation
        y <- as.numeric(res$Y[,k])
        rho <- cor( t(res$X), y)[,1]
        posf <- playbase::mofa.prefix(res$posf)

        par(mfrow=c(2,2), mar=c(4,4,2.5,1))
        if(length(res$posf)>4) par(mfrow=c(3,3))        
        for(i in 1:length(res$posf)) {
          pos1 <- posf[[i]]
          rho1 <- rho[rownames(pos1)]
          col1 <- playbase::colorscale(rho1, gamma=1)
          plot( pos1, col=col1, pch=20, cex=1.2, xlab="UMAP1", ylab="UMAP2" )
          title(toupper(names(res$posf)[i]), cex.main=1.2)          
        }
      }
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(90, 120),
      add.watermark = watermark
    )

    
  })
}
