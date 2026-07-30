##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Was: board.upload/upload_module_batchcorrect.R
if(0) {
  pgx <- playdata::GEIGER_PGX
}

qsee_imputation_server <- function(id, rX, rY) {
  shiny::moduleServer(
    id,
    function(input, output, session) {

      observeEvent( rY(), {
        cols <- colnames(rY())
        updateSelectInput(session, "colorby", choices = cols)
      })

      val_set <- NULL
      val_set <- reactiveVal(NULL)
      
      tr_marlevel <- reactive(input$marlevel) %>% debounce(1000)
      tr_mnarlevel <- reactive(input$mnarlevel) %>% debounce(1000)        
        
      get_normX <- reactive({
        message("[qsee_imputation_server] computing normX ")
        rawX <- rX()

        norm.method = "maxMedian"
        norm.method = "CPM+quantile"        
        normX <- playbase::normalizeExpression(rawX, method=norm.method, prior=0)        
        normX

        ##val_set <- NULL
        val_set(NULL)                

        ## add MNAR missing values
        mnarlevel = tr_mnarlevel()
        if(mnarlevel>0) {
          q001 <- quantile(normX, probs=0.01, na.rm=TRUE)
          sdx <- mean(matrixStats::colSds(scale(t(normX), scale=FALSE), na.rm=TRUE))          
          jj <- which(!is.na(normX))
          probx <- pnorm( -abs(normX[jj]-q001), mean=0, sd=sdx)
          ii <- sample(jj, mnarlevel * length(normX), prob=probx, replace=TRUE)
          ##val_set <- union(val_set, ii)
          vv <- cbind(ii, normX[ii])
          ##val_set <- rbind(val_set, vv) 
          val_set( rbind(isolate(val_set()), vv) )                    
          normX[ii] <- NA
        }

        ## add MAR missing values
        marlevel = tr_marlevel()
        if(marlevel>0) {
          jj <- sample( which(!is.na(normX)), marlevel * length(normX))
          vv <- cbind(jj, normX[jj])
          ##val_set <- rbind(val_set, vv) 
          val_set( rbind(isolate(val_set()), vv) )          
          normX[jj] <- NA
        }
        
        message(paste0("[qsee_imputation_server] missing values: ",
          round(100*mean(is.na(normX)),2),"%"))
        normX
      })

      get_impX <- reactive({
        X <- get_normX()        
        impX <- list("no imputation" = X)
        methods = c("SVD2", "QRILC", "MinProb", "Perseus")
        is.mox <- playbase::is.multiomics(rownames(X))
        m = "bpca"
        features <- shiny::withProgress(
          message = "Computing imputation methods...",
          for (m in methods) {
            if (is.mox) {
              impX[[m]] <- playbase::imputeMissing.mox(X, m)
            } else {
              impX[[m]] <- playbase::imputeMissing(X, m)
            }
          }
        )
        impX
      })

      get_pcaX <- reactive({
        impX <- get_impX()        
        pcaX <- list()
        features <- shiny::withProgress(
          message = "Computing PCA...",
          for (m in names(impX)) {
            cX <- impX[[m]] - rowMeans(impX[[m]], na.rm=TRUE)
            sel <- complete.cases(cX)
            res <- irlba::irlba(cX[sel,,drop=FALSE], nv=2)
            pcaX[[m]] <- res$v
            dimnames(pcaX[[m]]) <- list(colnames(cX),c("PC1","PC2"))
          }
        )
        pcaX
      })
      
      output$pca_plots <- renderPlot({        
        pcaX <- get_pcaX()

        ph <- input$colorby
        shiny::req(ph)
        y <- factor(rY()[,ph])
        cc <- 1 + as.integer(y)

        par(mfrow = c(2, 3), cex=1, mar = c(5, 5, 2, 1) )
        i=1
        for (m in names(pcaX)) {
          pos <- pcaX[[m]]
          plot(pos, main = m, pch=19, col = cc,
            xlab="PC2", ylab="PC2", cex=2)
          text(pos, labels = rownames(pos), cex=0.8, pos=3)
        }
      })


      output$heatmap_plots <- renderPlot({        
        impX <- get_impX()        

        ph <- input$colorby
        shiny::req(ph)        
        y <- factor(rY()[,ph])
        cc <- 1 + as.integer(y)
        
        cor.hclust <- function(x) {
          corx <- stats::cor(x, use = "pairwise")
          corx[is.na(corx)] <- 0
          fastcluster::hclust(stats::as.dist(1 - corx), method = "ward.D2")
        }

        names(impX)
        X1 <- impX[[1]]
        nax <- rowMeans(is.na(X1))
        sel <- ( nax < 0.50 & nax > 0)
        X1 <- X1[sel,,drop=FALSE]
        X1 <- X1[order(-matrixStats::rowSds(X1,na.rm=TRUE)),]
        X1 <- head(X1,200)
        gg <- rownames(X1)
        ii <- cor.hclust(t(X1))$order
        jj <- cor.hclust(X1)$order
        
        par(mfrow = c(2, 3), cex=1, mar = c(5, 1, 2, 6) )
        i=1
        for (m in names(impX)) {
          X1 <- impX[[m]][gg,]
          X1 <- X1[ii,jj]
          playbase::gx.imagemap(X1, clust=FALSE, main=m, cex.main=1.2)
        }

      })

      
      output$histograms <- renderPlot({
        impX <- get_impX()        
        
        ii <- which(!is.na(impX[[1]]))
        jj <- which(is.na(impX[[1]]))        
        
        par(mfrow = c(2, 3), cex=1, mar = c(5, 5, 3, 1) )
        m=names(impX)[2]
        for (m in names(impX)) {
          X1 <- impX[[m]]
          hist( X1[ii], breaks=80, main=m, border="grey60",
            xlab="intensity (log2)")
          if(any(!is.na(X1[jj]))) {
            hist( X1[jj], breaks=100, add=TRUE, col='red', border=NA)
          }
        }
        
      })

      output$distributions <- renderPlot({
        impX <- get_impX()        
        
        ii <- which(!is.na(impX[[1]]))
        jj <- which(is.na(impX[[1]]))        
        
        ph <- input$colorby
        shiny::req(ph)
        y <- factor(rY()[,ph])
        cc <- 1 + as.integer(y)

        par(mfrow = c(2, 3), cex=1, mar = c(5, 5, 3, 1) )

        X1 <- impX[[1]]
        x.avg <- rowMeans(X1, na.rm = TRUE)
        x.nar <- rowMeans(is.na(X1))
        x.avg2 <- cut(x.avg, breaks = 20)
        x.nar2 <- tapply(1:nrow(X1), x.avg2, function(i) mean(is.na(X1[i, , drop = FALSE])))
        aa <- sort(unique(as.numeric(gsub(".*,|\\]", "", as.character(x.avg2)))))
        R <- rbind(x.nar2, 1 - x.nar2)
        barplot(
          R,
          beside = FALSE, names.arg = aa, las = 1,
          col = c("red","lightgrey"),
          border = "grey60", 
          xlab = "average intensity (log2)",
          ylab = "missing ratio"
        )
        title("missingness vs. intensity", cex.main=1.2)

        plot(x.avg, x.nar,
          xlab = "average intensity (log2)",
          ylab = "missing ratio")          
        title("missingness vs. intensity", cex.main=1.2)
        
        hist( x.nar[x.nar>0], breaks=40,
          main = "missing ratio histogram",
          xlab = "missing ratio", cex.main=1.2)

        qq <- seq(0,1,0.05)
        qsum <- sapply(qq, function(q) sum(x.nar<=q))
        plot(qq, qsum,
          xlab = "max. missing ratio threshold",
          ylab = "nr. features")
        title("nr. features vs. threshold", cex.main=1.2)
        
        par(mar = c(7, 5, 3, 1) )
        barplot( colMeans(is.na(X1)),
          las = 3, cex.names = 0.85, col = cc,
          ylab = "missing ratio")          
        title("missingness per sample", cex.main=1.2)        
        
      })

      output$validation_scatter <- renderPlot({
        impX <- get_impX()        
        V <- val_set()
        if(length(V)==0) {
          par(mfrow=c(1,1))
          plot.new();
          text(0.5,0.5,"Please add MAR/MNAR missing values.")
          return(NULL)
        }

        i=1
        par(mfrow=c(2,3), cex=1.1, mar=c(4,4,3,1), mgp=c(2.5,0.8,0))
        for(i in 2:length(impX)) {
          X1 <- impX[[i]]
          jj <- V[,1]
          plot( X1[jj], V[,2],
            pch=19, cex=0.4,
            xlim=c(-5,15),
            main = names(impX)[i],
            xlab = "imputed value",
            ylab = "actual value"
          )
          rho <- cor(X1[jj], V[,2])
          tt <- paste("R =", round(rho,2))
          legend("bottomright", tt, bty='o', bg="white", cex=0.9)
        }
        
        
      })
      
      ## Keep this output alive even when the Normalization tab is hidden
      ##outputOptions(output, "normalize_plots", suspendWhenHidden = FALSE)
      
    } ## end-of-server
  )
}
