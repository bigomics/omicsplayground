##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_lasagna_clustering_ui <- function( id,  ... )
{
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::radioButtons(ns("colorby"), "Color by:",
                        c("foldchange","correlation"),
                        selected="correlation")
  )
  
  PlotModuleUI(
    ns("plot"),
    download.fmt = c("png", "pdf", "svg"),
    options = options,
    ...
  )
}

mofa_plot_lasagna_clustering_server <- function(id,
                                                data,
                                                input_contrast = reactive(NULL),
                                                watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {
      res <- data()
      k <- input_contrast()
      shiny::req(res$posx)
      shiny::req(res$posf)
      shiny::req(res$Y)      
      if(!is.null(k)) shiny::req(k %in% colnames(res$Y))
      
      col1 = 'black'      
      y <- as.numeric(res$Y[,k])
      y <- sign(y)
      y[y==0] <- NA
      if(input$colorby=="correlation") {
        ## color by contrast correlation
        rho <- cor(t(res$X), y, use="pairwise")[,1]
      } else {
        ## color by fold-change
        m1 <- rowMeans(res$X[, which(y == 1)],na.rm=TRUE)
        m0 <- rowMeans(res$X[, which(y == -1)],na.rm=TRUE)
        rho <- m1 - m0
      }
      posf <- playbase::mofa.prefix(res$posf)
      posf <- posf[setdiff(names(posf),c("PHENO"))]
      
      par(mfrow=c(2,2), mar=c(4,4,2.5,1))
      if(length(posf)>4) par(mfrow=c(3,3))        
      for(i in 1:length(posf)) {
        pos1 <- posf[[i]]
        rho1 <- rho[rownames(pos1)]
        col1 <- playbase::colorscale(rho1, gamma=1)
        xlab <- colnames(pos1)[1]
        ylab <- colnames(pos1)[2]
        plot(pos1, col = col1, pch = 20, cex = 1.4,
          xlab = xlab, ylab = ylab, las = 1)
        title(toupper(names(posf)[i]), cex.main = 1.2)          
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
