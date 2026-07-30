##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Was: board.upload/upload_module_batchcorrect.R
if(0) {
  pgx <- playdata::GEIGER_PGX
}

qsee_pca_server <- function(id, rX, rY) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      
      get_normX <- reactive({
        message("[normalize_plots] computing normX ")
        rawX <- rX()
        m = "CPM+quantile"
        normX <- playbase::normalizeExpression(rawX, method = m, ref=NULL, prior=0)
        return(normX)
      })

      get_pcaX <- reactive({
        normX <- get_normX()        
        cX <- normX - rowMeans(normX, na.rm=TRUE)
        sel <- which(rowMeans(is.na(cX))==0)  ## only complete rows
        nv <- min(10, min(dim(cX))-1)
        pca <- irlba::irlba(cX[sel,], nv=nv)
        rownames(pca$u) <- rownames(cX)[sel]
        rownames(pca$v) <- colnames(cX)

        Y <- rY()
        H <- playbase::expandPhenoMatrix(Y, drop.ref=FALSE)
        pca$R <- cor(H, pca$v)

        pca
      })
      
      output$pairs <- renderPlot({        
        pca <- get_pcaX()        
        y <- rY()[,1]
        cc <- 1 + as.integer(factor(y))
        V <- pca$v[,1:4]
        colnames(V) <- paste0("PC",1:ncol(V))
        pairs(V, col=cc, pch=19, cex=2,
          oma=c(12,5,6,1), main="PC scatterplot matrix", line.main=3)
      })

      output$variance <- renderPlot({        
        pca <- get_pcaX()
        aa <- paste0("PC",1:ncol(pca$v))
        dd <- 100 * pca$d**2 / sum(pca$d**2)
               
        par(cex=1, mar=c(4.5,4,3,1))
        barplot(dd, names.arg=aa,
          ylab = "variance (%)", width=1,
          xlab = "PC component",
          ylim = c(0,max(dd)*1.1)          
        )
        title("Variance explained")
        xx <- (1:length(dd) - 0.5) * 1.2
        lines(xx, dd, type="b", cex=1.6, pch=20)
        text(xx, dd, labels=paste0(round(dd,2),"%"), pos=3, cex=1)
        if(0) {
          cdd <- cumsum(dd)
          lines(xx, cdd, type="b", cex=1.6, pch=20)        
          text(xx, cdd, labels=paste0(round(cdd,2),"%"), pos=3, cex=1)
        }
      })

      output$trait_correlation <- renderPlot({        
        pca <- get_pcaX()
        Y <- rY()
        H <- playbase::expandPhenoMatrix(Y, drop.ref=FALSE)
        V <- pca$v
        colnames(V) <- paste0("PC",1:ncol(V))
        R <- cor(V, H)

        px <- sub("=.*","",colnames(R))
        B <- tapply(1:ncol(R), px, function(i) rowMeans(R[,i,drop=FALSE]**2))
        B <- do.call(cbind, B)
        B <- B / rowSums(B)
        
        ##heatmap(R)
        par(cex=1, mar=c(4,1,3,10))
        playbase::gx.imagemap(t(R), clust=FALSE, col=heat.colors(32),
          main="PC-trait correlation", cex.main=1.2)                
      })
      
      output$feature_pca <- renderPlot({        
        pca <- get_pcaX()
        normX <- get_normX()
        
        pos <- pca$u[,1:2]
        ## pos <- uwot::umap2( pca$u )        
        Y <- rY()
        H <- playbase::expandPhenoMatrix(Y, drop.ref=FALSE)        
        rx <- cor(t(normX[rownames(pos),]), H, use="pairwise")

        par(mfrow=c(3,3), mar=c(3.8,4,2.2,1), mgp=c(2.4,0.8,0))
        if(ncol(rx)>9) par(mfrow=c(4,4))
        for(i in 1:min(ncol(rx),16)) {
          cc <- playbase::colorscale(rx[,i])
          plot( pos, col=cc, #main = colnames(rx)[i],
            xlab="PC1", ylab="PC2", pch=19 )
          title(main = colnames(rx)[i], line=0.8)
          
          dir <- colMeans(pos * pmax(rx[,i],0, na.rm=TRUE)**2 )
          dir <- 0.05 * dir / sqrt(sum(dir**2))
          arrows( 0, 0, dir[1], dir[2], length=0.15, lwd=3 )
        }

      })
      
      output$biplot <- renderPlot({        
        res <- get_pcaX()
        samples <- rY()
        par(cex=1.1, mar=c(12,4,5,0.8))
        arrow_type <- input$arrows
        maxarrows <- as.integer(input$maxarrows)
        pca.draw_biplot(res, samples, xpc=1, ypc=2,
          main = "PCA biplot", cex = 1.2, cex.labels = 0.9,          
          colorby_var=1, arrows=arrow_type, maxarrows=maxarrows) 
      })

      output$biplot2 <- renderPlot({
        require(ggplot2)
        
        X <- get_normX()
        Y <- rY()
        groups <- Y[,1]
        H <- playbase::expandPhenoMatrix(Y, drop.ref=FALSE)        
        
        X1 <- X
        #X1 <- playbase::svdImpute2(X1, nv=5)
        sel <- which(rowSums(is.na(X1))==0)  ## only complete rows??
        xx <- rbind( t(H)*0.001, X1[sel,])
        res <- prcomp( t(xx), scale. = FALSE)

        type = "loadings"
        type = "pheno"
        type <- input$arrows
        if(type == "pheno") {
          U <- res$rotation[1:ncol(H),,drop=FALSE]
        } else {
          U <- res$rotation[-(1:ncol(H)),,drop=FALSE]        
        }
        
        dims <- c(1,2)
        U <- scale(U, center=FALSE)        
        U <- U[order(-rowSums(U[,dims]**2)),]
        
        maxarrows=10
        maxarrows <- as.integer(input$maxarrows)
        U <- head(U, maxarrows)
        res$rotation <- U
        #res$sdev <- rep(1,ncol(U))
        #res$x <- scale(res$x, center=FALSE)
        
        cex = 1.2
        par(mar=c(4,4,2,2))
        g <- ggbiplot::ggbiplot(
          res,
          choices = dims,
          point.size = 2.2*cex,
          obs.scale = 1,
          var.scale = 1, 
          varname.size = 3.8*cex,
          labels = rownames(Y),
          labels.size = 3.5*cex,
          groups = groups, 
          ellipse = TRUE,
          circle = TRUE
        ) +
          ggtitle("PCA biplot") +
          ggplot2::labs(fill = "Groups", color = "Groups") +
          ggplot2::geom_point(ggplot2::aes(color = groups))
        
        
        g <- g +
          ggplot2::theme(
            axis.text = ggplot2::element_text(size=10*cex),
            axis.title = ggplot2::element_text(size=12*cex),
            plot.margin = unit(c(1,1,1,1)*0.4, "cm")
          )

        #g + ggprism::theme_prism(base_size = 18)
        g <- g + theme_light(base_size = 12*cex)
        #g <- g + theme_bw(base_size = 12*cex)
        #g + theme_grey(base_size = 18) 

        g + ggplot2::theme(
          aspect.ratio = 1,
          legend.direction = 'horizontal', legend.position = 'bottom',
        )
        
      })
      
    } ## end-of-server
  )
}


## Biplot grid for the BC panel: one biplot per method (sample scores as
## points, arrows overlaid) on the two selected PCs. Same method list /
## grid as plot_all_methods. Two arrow sources:
##  - "loadings": top feature loadings (res$loadings, same PCA as scores)
##  - "pheno": each annotation's correlation with the two PCs, precomputed
##    in playbase (res$pheno.cor; eigencorplot-style, plain cor())
pca.draw_biplot <- function(res, samples, xpc=1, ypc=2, colorby_var=NULL,
                            main = "PCA biplot", cex = 2, cex.labels = 0.7,
                            arrow.col = "#7A2E2E", arrow.lwd = 2,
                            arrows = "loadings", maxarrows=10) {

  if(0) {
    X <- playbase::logCPM(playbase::COUNTS)
    X <- X - rowMeans(X,na.rm=TRUE)
    Y <- playbase::SAMPLES
    res <- svd(X)
    rownames(res$v) <- rownames(Y)
    rownames(res$u) <- rownames(X)
    H <- playbase::expandPhenoMatrix(Y, drop.ref=FALSE)
    res$R <- cor(H, res$v)
    xpc=1;ypc=2
  }
  
  ## arrow endpoints for a method: named 2-col matrix in loading- or
  ## correlation-units (scaled to the score cloud later). Both come
  ## precomputed from playbase::compare_batchcorrection_methods.
  arrow_ends <- function(res, scores, arrows="loadings", maxarrows=10) {
    if (arrows == "loadings") {
      L <- res$u
      if (is.null(L) || max(xpc, ypc) > ncol(L)) {
        return(NULL)
      }
      A <- L[, c(xpc, ypc), drop = FALSE]
    } else {
      Rc <- res$R
      if (is.null(Rc) || max(xpc, ypc) > ncol(Rc)) {
        return(NULL)
      }
      A <- Rc[, c(xpc, ypc), drop = FALSE]
      if (!nrow(A)) {
        return(NULL)
      }
    }
    A <- A[stats::complete.cases(A), , drop = FALSE]      
    A[head(order(-(A[, 1]^2 + A[, 2]^2)), maxarrows), , drop = FALSE]
  }
  
  P <- res$v
  if (is.null(P) || max(xpc, ypc) > ncol(P)) {
    plot.new()
    text(0.45, 0.5, "method failed")
    title(m, cex.main = 1.3)
    return(invisible())
  }
  scores <- P[, c(xpc, ypc), drop = FALSE]
  smp <- samples[rownames(scores), , drop = FALSE]
  color <- 1 + as.integer(factor(smp[, colorby_var]))

  cex1 <- cut(nrow(scores), c(0, 40, 100, 250, 1000, 999999),
    c(1, 0.85, 0.7, 0.55, 0.4))
  cex1 <- 2 * cex * as.numeric(as.character(cex1))
  if (is.na(cex1)) cex1 <- 1

  pcvar <- res$d**2 / sum(res$d**2)
  pcvar <- round( 100 * pcvar[c(xpc, ypc)], 1)
  
  lim <- c(-1, 1) * max(abs(scores)) * 1.15
  plot( scores,
    col = color, pch = 20, cex = cex1, las = 1,
    xlim = lim, ylim = lim,
    xlab = paste0("PC", xpc, " (", pcvar[1], "%)"),
    ylab = paste0("PC", ypc, " (", pcvar[2], "%)"),
  )
  title(main, cex.main = 1.3, line=1)
  abline(h = 0, v = 0, lty = 3, col = "grey70")

  plotrix::thigmophobe.labels(scores[,1], scores[,2], rownames(scores),
    cex=cex.labels, col='black', offset = 0.4, xpd = NA
  )
  
  A <- arrow_ends(res, scores, arrows = arrows, maxarrows = maxarrows)
  if (is.null(A) || !nrow(A) || max(abs(A)) == 0) {
    return(invisible())
  }
  ## scale arrows to fill the score cloud (loadings/correlations sit on a
  ## much smaller scale than the scores)
  s <- 0.85 * max(abs(scores)) / max(abs(A))
  ax <- A[, 1] * s
  ay <- A[, 2] * s
  arrows(0, 0, ax, ay, length = 0.15, col = arrow.col, lwd = arrow.lwd)
  ## repel labels off each other (base-graphics; places each away from
  ## its nearest neighbour) rather than a fixed left/right offset
  plotrix::thigmophobe.labels(ax, ay, rownames(A),
    cex = cex.labels, col = arrow.col, offset = 0.4, xpd = NA
  )
}


