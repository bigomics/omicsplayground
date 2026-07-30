##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Was: board.upload/upload_module_batchcorrect.R
if(0) {
  pgx <- playdata::GEIGER_PGX
}

qsee_filtering_server <- function(id, rX, rY) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      
      shiny::observeEvent(rY(), {
        cols <- colnames(rY())
        shiny::updateSelectInput(
          session, "colorby",
          choices = cols,
          selected = cols[1]
        )
      }, ignoreNULL = TRUE)

      shiny::observeEvent(rX(), {
        X <- rX()
        nmax <- max(1L, nrow(X))
        n0 <- min(1000L, nmax)
        shiny::updateSliderInput(
          session, "threshold",
          min = 1,
          max = nmax,
          value = n0,
          step = 1
        )
      }, ignoreNULL = TRUE)

      get_pcaX <- shiny::reactive({
        X <- rX()
        Y <- rY()

        ## only complete rows. REALLY???
        sel <- which(rowSums(is.na(X)) == 0)  
        X <- X[sel,]
        
        cX <- X - rowMeans(X, na.rm = TRUE)  ## important
        nv <- min(10, dim(cX) - 1)
        res <- irlba::irlba(cX, nv = nv)
        rownames(res$u) <- rownames(cX)
        rownames(res$v) <- colnames(cX)
        colnames(res$u) <- paste0("PC", 1:ncol(res$u))
        colnames(res$v) <- paste0("PC", 1:ncol(res$v))

        H <- playbase::expandPhenoMatrix(Y, drop.ref = FALSE)
        res$R <- stats::cor(H, res$v)

        res$X <- X
        res$Y <- Y
        res$H <- H

        res
      })
      
      output$pca_vs_topsd <- shiny::renderPlot({
        res <- get_pcaX()
        ph <- input$colorby
        shiny::req(res, ph)
        pcaexplorer.plot_pca_vs_topsd(res, colorby_var = ph)
      })

      output$variance_vs_topsd <- shiny::renderPlot({
        res <- get_pcaX()
        top <- input$threshold
        shiny::req(res, top)
        pcaexplorer.plot_variance_vs_topsd(res, topsd = top)
      })

      output$sd_histogram <- shiny::renderPlot({
        res <- get_pcaX()
        top <- input$threshold
        shiny::req(res, top)
        pcaexplorer.plot_sd_histogram(res, topsd = top)
      })


      
    } ## end-of-server
  )
}
