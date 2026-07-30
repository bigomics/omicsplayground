##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Was: board.upload/upload_module_batchcorrect.R
if(0) {
  pgx <- playdata::GEIGER_PGX
}

qsee_normalization_server <- function(id, rX, rY) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      
      observeEvent( rY(), {
        cols <- colnames(rY())
        updateSelectInput(session, "colorby", choices = cols)
      })

      ## Inject technical noise; amount is controlled by the sidebar slider
      ## (0 = none, 1 = default level used previously in qsee_server).
      get_rawX <- reactive({
        rawX <- rX()
        shiny::req(rawX)
        amount <- input$noise
        if (is.null(amount)) amount <- 1
        if (amount <= 0) {
          return(rawX)
        }
        ii <- which(rawX == 0)
        n <- ncol(rawX)
        ## Stable draws for a given noise level so plots don't jump on redraw
        set.seed(42 + as.integer(round(amount * 10)))
        ## Multiplicative + additive technical variation, scaled by amount
        rawX <- t(t(rawX) * rnorm(n, mean = 1, sd = 0.2 * amount) +
          10 * amount + rnorm(n, mean = 0, sd = 2 * amount))
        rawX <- pmax(rawX, 0)
        rawX[ii] <- NA
        rawX
      })

      get_normX <- reactive({
        rawX <- get_rawX()
        ## methods (matching those used in the full upload module)
        methods <- c("CPM", "CPM+quantile", "TMM", "maxMedian", "maxSum", "quantile")
        methods <- c("CPM", "CPM+quantile", "maxMedian", "maxSum", "quantile")                
        normX <- list(raw = rawX)
        prior <- 0
        for (m in methods) {
          normX[[m]] <- tryCatch(
            playbase::normalizeExpression(rawX, method = m, ref = NULL, prior = prior),
            error = function(e) {
              message("normalizeExpression failed for ", m, ": ", e$message)
              rawX  ## fallback to raw data
            }
          )
        }
        return(normX)
      })

      get_pcaX <- reactive({
        normX <- get_normX()        
        pcaX <- list()
        for (m in names(normX)) {
          cX <- normX[[m]] - rowMeans(normX[[m]], na.rm=TRUE)
          sel <- which(rowMeans(is.na(cX))==0)  ## only complete rows
          res <- irlba::irlba(cX[sel,], nv=2)
          pcaX[[m]] <- res$v
          dimnames(pcaX[[m]]) <- list(colnames(cX),c("PC1","PC2"))
        }        
        pcaX
      })
      
      output$box_plots <- renderPlot({        
        normX <- get_normX()
        ph <- input$colorby
        shiny::req(normX, ph)

        y <- factor(rY()[,ph])
        cc <- 1 + as.integer(y)
        par(mfrow = c(2, 3), cex=1, mar = c(7, 5, 3, 1), las=2)
        for (i in seq_along(normX)) {
          x <- normX[[i]]
          if (is.null(x) || length(x) == 0) {
            plot(0, type = "n", main = names(normX)[i])
            text(0, 0, "No data")
            next
          }
          ## Sample columns if too many samples (common in real data)
          if (ncol(x) > 40) {
            x <- x[, sample(ncol(x), 40)]
          }
          boxplot(x, col = cc, main = names(normX)[i], ylab = "log-expression")
        }
      })

      output$histograms <- renderPlot({        
        normX <- get_normX()
        ph <- input$colorby
        shiny::req(normX, ph)

        y <- factor(rY()[,ph])
        cc <- 1 + as.integer(y)

        par(mfrow = c(2, 3), cex=1, mar = c(6, 5, 3, 1), las=2)
        i=1
        for (i in seq_along(normX)) {
          x <- normX[[i]]
          if (is.null(x) || length(x) == 0) {
            plot(0, type = "n", main = names(normX)[i])
            text(0, 0, "No data")
            next
          }
          ## Sample columns if too many samples (common in real data)
          if (ncol(x) > 40) {
            x <- x[, sample(ncol(x), 40)]
          }
          playbase::gx.hist(x, main=names(normX)[i], breaks=80,
            lines.col=cc, las=1, xlab = "signal")
        }

      })
      
      output$pca_plots <- renderPlot({        
        pcaX <- get_pcaX()
        ph <- input$colorby
        shiny::req(pcaX, ph)
        
        y <- factor(rY()[,ph])
        cc <- 1 + as.integer(y)
        par(mfrow = c(2, 3), cex=1, mar = c(5, 5, 3, 1) )
        for (i in seq_along(pcaX)) {
          pos <- pcaX[[i]]
          plot(pos, main = names(pcaX)[i], pch=19, col = cc,
            xlab="PC1", ylab="PC2", cex=2)
          text(pos, labels = rownames(pos), cex=0.8, pos=3)
        }
      })
      
    } ## end-of-server
  )
}
