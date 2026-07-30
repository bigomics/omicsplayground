##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Was: board.upload/upload_module_batchcorrect.R
if(0) {
  pgx <- playdata::GEIGER_PGX
}

qsee_outlier_server <- function(id, rX, rY) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      
      observeEvent( rY(), {
        cols <- colnames(rY())
        updateSelectInput(session, "colorby", choices = cols)
      })

      output$outlier_zscores <- renderPlot({
        X <- rX()
        ph <- input$colorby
        shiny::req(X, ph)

        X <- X[complete.cases(X),,drop=FALSE]
        y <- rY()[,ph]
        cc <- 1 + as.integer(factor(y))

        par(mfrow = c(2,3), mar = c(8, 4.5, 2, 1), cex=1, mgp=c(2.3,0.8,0) )        
        methods <- c("z.correlation", "z.distance", "z.features", "z.isoforest")
        out <- playbase::detectOutlierSamples(X, methods=methods,
          col=cc, plot=TRUE, par=FALSE)
      })

      output$outlier_heatmap <- renderPlot({
        X <- rX()
        shiny::req(X)
        playbase::gx.heatmap(X, nmax=400, keysize=1, mar=c(12,8) )
      })

      output$outlier_pca <- renderPlot({
        X <- rX()
        ph <- input$colorby
        shiny::req(X,ph)

        y <- rY()[,ph]
        cc <- 1 + as.integer(factor(y))

        cX <- X - rowMeans(X,na.rm=TRUE)
        if(any(is.na(cX))) cX <- playbase::svdImpute2(cX)
        res <- svd(cX)
        plot( res$v[,1:2], col=cc, pch=19, cex=2.4,
          xlab = "PC1", ylab = "PC2")
        
      })
      
    } ## end-of-server
  )
}
