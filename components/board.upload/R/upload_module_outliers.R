##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## =============================================================================
## ==================== OUTLIERS UI/SERVER =================================
## =============================================================================


upload_module_outliers_ui <- function(id, height = "100%") {
  ns <- shiny::NS(id)
  
  score.infotext =
  "The outlier z-score is calculated as the average of z-score from correlation, euclidean distance and avarage feature z-score."

  missing.infotext =
  "Analysis of variables by plotting their significance in correlation with the phenotype against their significance in correlation with a principal component (PC) vector. Strong model variables are situate 'top right'. Batch effect variables with high PC correlation but low phenotype correlation are on the 'top left'. A well-designed experiment shows strong model variables in PC1, else it may be a sign of significant batch-effects."
  
  norm.options <- tagList(
    shiny::radioButtons( ns('norm_plottype'), "Plot type:", c("boxplot"),
      selected = "boxplot", inline = TRUE),
  )

  outlier.options <- tagList(
    shiny::radioButtons( ns('outlier_plottype'), "Plot type:", c("pca","heatmap"), inline = TRUE),
    shiny::checkboxInput( ns('outlier_shownames'), "show sample names", FALSE)
  )

  bec.options <- tagList(
    shiny::radioButtons( ns('bec_plottype'), "Plot type:", c("pca","tsne","umap","heatmap"),
      inline = TRUE)
  )

  bslib::layout_columns(
    col_widths = c(2,10),
    height = "calc(100vh - 200px)",
    heights_equal = "row",
    bslib::card( bslib::card_body(
      style = "padding: 0px;",
      bslib::accordion(
        multiple = FALSE,
        style = 'background-color: #F7FAFD99;',
        bslib::accordion_panel(      
          title = "1. Normalization",
          div("Normalize data values:\n"),
          shiny::selectInput( ns("scaling_method"), NULL,
            choices = c("CPM"="cpm","median.4"="m4","sd.001" = "sd"),
            selected = "cpm"
          ),
          shiny::checkboxInput(ns("clip_zero"), "clip to zero", value=TRUE),
          shiny::checkboxGroupInput(ns("correct_options"), "Correct for:", ## inline = TRUE,
            choices = c(
              "library effects" = "lib",
              "cell state (mito/ribo/CC)" = "cell",
              "gender (M/F)" = "gender"
            ),
            selected = c("lib","cell","gender")),
          br()          
        ),      
        bslib::accordion_panel(      
          title = "2. Missing values",
          shiny::p("Replace missing values using an imputation method:\n"),
          shiny::selectInput( ns("impute_method"), NULL,
            ##  choices = c("bpca","LLS","MinDet","MinProb","NMF","RF","SVD2","zero"),
            choices = c("MinDet","MinProb","NMF","SVD2","zero"),
            selected='SVD2'),
          shiny::checkboxInput( ns("zero_as_na"), label="Treat zero as NA", value=TRUE),
          br()          
        ),
        bslib::accordion_panel(
          title = "3. Outlier detection",
          p("Analyze and remove outliers (i.e. bad samples) from your dataset.\n"),
          shiny::sliderInput( ns("outlier_threshold"), "Threshold:", 1, 6, 6, 1),
          br()          
        ),
        bslib::accordion_panel(      
          title = "4. Batch correction",
          shiny::p("Clean up your data from unwanted variation:\n"),
          shiny::selectInput( ns("bec_method"), NULL,
            choices = c("<none>","ComBat","limma","SVA","RUV3","NNM"),
            selected = "SVA"
          ),
          shiny::conditionalPanel(
            "input.bec_method == 'ComBat' || input.bec_method == 'limma'",
            ns = ns,
            shiny::selectInput(ns("bec_param"), "Batch parameter:", choices=NULL)            
          ),
          shiny::checkboxInput(ns("bec_preview"), "Preview all methods", value = FALSE),          
          shiny::actionButton(
            ns("compute_button"),
            "Compute",
            class = "btn btn-primary",
            width ='80%',
            style = 'margin-left:10%; margin-top: 15px;'
            ##style = 'margin-left:10%; position: absolute; bottom: 30px;'
          ),
          br()
        )
      ),
      br()
    )),

    ## ---------------------------- canvas ----------------------------------
    bslib::layout_columns(
      width = 12,
      bslib::layout_columns(
        col_widths = 6,
        row_heights = c(3,3),
        height = "calc(100vh - 200px)",
        heights_equal = "row",
        ##  shiny::plotOutput(ns("canvas"), width = "100%", height = height) %>% bigLoaders::useSpinner(),
        PlotModuleUI(
          ns("plot1"),
          title = "Normalization",
#         info.text = info.text,
#         caption = caption,
          options = norm.options,
          height = c("100%","70vh")          
        ),        
        PlotModuleUI(
          ns("plot2"),
          title = "Missing values",
          info.text = missing.infotext,
          caption = missing.infotext,
          options = NULL,
          height = c("100%","70vh")
        ),
        PlotModuleUI(
          ns("plot3"),
          title = "Outlier detection",
          info.text = score.infotext,
          caption = score.infotext,
          options = outlier.options,
          height = c("100%","70vh")          
        ),
        PlotModuleUI(
          ns("plot4"),
          title = "Batch-effect correction",
#          info.text = info.text,
#          caption = caption,
          options = bec.options,
          height = c("100%","70vh")           
        )
      )
    )
  )
}


upload_module_outliers_server <- function(id, r_X, r_samples, r_contrasts,
                                              is.count = FALSE, height = 720) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      

      observeEvent({
        list( get_batchpars() )
      }, {       
        pars <- get_model_parameters()
        batch.pars <- pars$batch
        shiny::req( batch.pars )
        shiny::updateSelectInput( session, "bec_param", choices = batch.pars)
        
      })
      
      get_model_parameters <- eventReactive({
        list( r_X(), r_samples(), r_contrasts() )
      }, {
        
        shiny::req( r_X(), r_samples(), r_contrasts() )

        X <- r_X()
        samples <- r_samples()
        contrasts <- r_contrasts()
        ## contrasts <- contrasts[,1,drop=FALSE]
        
        dbg("[outliers_server:observeEvent] 1: dim.X = ", dim(X) )
        dbg("[outliers_server:observeEvent] 1: dim.samples = ", dim(samples) )
        dbg("[outliers_server:observeEvent] 1: dim.contrasts = ", dim(contrasts) )

        shiny::req( ncol(X) == nrow(samples) )
        shiny::req( nrow(contrasts) == ncol(X) )

        dbg("[outliers_server:observeEvent] 2: dim.X = ", dim(X) )
        dbg("[outliers_server:observeEvent] 2: dim.samples = ", dim(samples) )
        dbg("[outliers_server:observeEvent] 2: dim.contrasts = ", dim(contrasts) )
        
        pheno <- playbase::contrasts2pheno( contrasts, samples )        
        bc <- playbase::detectBatchEffects(X, samples, pheno,  params = "statistical",
          k.pca=10, p.pca=0.5, p.pheno=0.05, xrank = 10)

        pheno.pars <- names(which(bc$p.values[,'p.pheno'] < 1e-8))
        pheno.pars

        batch.pars <- NULL
        if("statistical" %in% names(bc$params)) {
          batch.pars <- bc$params$statistical
          batch.pars
          batch.pars <- setdiff(batch.pars, pheno.pars)
          batch.pars2 <- grep("batch",colnames(samples),ignore.case=TRUE,value=TRUE)
          if(length(batch.pars2)) batch.pars <- c(batch.pars, batch.pars2)
          batch.pars <- sort(batch.pars)
        }

        batch.vec <- NULL
        if(length(batch.pars)) {
          batch.vec <- apply(samples[,batch.pars],1,paste,collapse='_')
        }
        pheno.vec <- NULL
        if(length(pheno.pars)) {
          pheno.vec <- apply(samples[,pheno.pars],1,paste,collapse='_')
        }
        
        list( batch.pars = batch.pars, pheno.pars = pheno.pars,
          batch.vec = batch.vec, pheno.vec = pheno.vec)
      })
      
      ##------------------------------------------------------------------
      ## Object reactive chain
      ##------------------------------------------------------------------

      normalizedX <- reactive({
        X <- r_X()
        samples <- r_samples()
        contrasts <- r_contrasts()
        shiny::req(X, samples)

        X <- playbase::pgx.countNormalization(X, methods = "median.center", keep.zero = TRUE)
        X <- playbase::scale_counts(X, method = input$scaling_method)
        logX <- log2(1 + X)

        ## technical effects correction
        opts <- input$correct_options
        dbg("[normalizedX] cell correction options = ",opts)        
        if(length(opts)) {
          pheno <- playbase::contrasts2pheno( contrasts, samples)          

          dbg("[normalizedX] head.pheno = ",head(pheno))

          bc <- playbase::detectBatchEffects( logX, samples, pheno,  params = "technical",
            k.pca=10, p.pca=0.5, p.pheno=0.05, xrank = 100)

          B <- bc$covariates
          sel <- c()
          if("lib" %in% opts) {
            sel <- c(sel, grep("^lib", colnames(B)))
          }
          if("cell" %in% opts) {
            sel <- c(sel, grep("^mito", colnames(B)))
            sel <- c(sel, grep("^ribo", colnames(B)))
            sel <- c(sel, grep("^cellcycle", colnames(B)))
          }
          if("gender" %in% opts) {
            sel <- c(sel, grep("^gender", colnames(B)))
          }
          B <- B[,sel,drop=FALSE]
          dbg("[normalizedX] dim(B) = ",dim(B))
          mod1 <- model.matrix(~pheno)
          B <- scale(B)
          B[is.na(B)] <- 0
          dbg("[normalizedX] *** correcting for  colnames(B) = ",colnames(B))
          dbg("[normalizedX] dim(B) = ",dim(B))
          res <- try(limma::removeBatchEffect(logX, covariates=B, design = mod1))
          if(!"try-error" %in% class(res)) {
            logX <- res
          }
        }

        ## for quantile normalization we omit the zero value and put back later
        jj <- which( logX < 0.01 )
        logX[jj] <- NA
        logX <- limma::normalizeQuantiles( logX )
        rownames(logX) <- rownames(X)
        colnames(logX) <- colnames(X)
        logX[jj] <- 0
        logX
      })

      imputedX <- reactive({
        shiny::req(normalizedX())
        offset <- 1e-20
        X <- normalizedX()
        if(input$clip_zero) X <- pmax(X, 0)
        X[playbase::is.xxl(X, z = 10)] <- NA  ## outlier XXL values
        if(input$zero_as_na)  X[which(X==0)] <- NA
        ##which.missing <- which( is.na(X) )
        X <- playbase::imputeMissing(X, method = input$impute_method )
        X <- pmax(X, 0)
        X
      })
      
      cleanX <- reactive({
        shiny::req(imputedX())        
        X <- imputedX()
        res <- playbase::detectOutlierSamples(X, plot=FALSE, y=NULL)        
        is.outlier <- ( res$z.outlier > input$outlier_threshold )
        if(any(is.outlier) && !all(is.outlier) ) {
          X <- X[, which(!is.outlier), drop = FALSE]
        }
        pos <- NULL
        if(NCOL(X)>1) {
          dbg("[cleanX:irlba] dim.X = ", dim(X))
          pos <- irlba::irlba(X, nv=2)$v
          rownames(pos) <- colnames(X)        
        }
        list( X=X, pos=pos )
      })
      
      correctedX <- shiny::eventReactive({
        input$compute_button
      },{
        ## triggers too often.... need button??
        shiny::req(cleanX()$X, r_contrasts(), r_samples())                
        X <- cleanX()$X
        
        samples <- r_samples()
        contrasts <- r_contrasts()
        batch.param <- input$bec_param
        
        samples <- samples[colnames(X),]
        contrasts <- contrasts[colnames(X),,drop=FALSE]        
        
        pheno <- playbase::contrasts2pheno( contrasts, samples )
        pheno <- pheno[colnames(X)]        
        m <- input$bec_method        

        dbg("[correctedX] length.batch.param = ", length(batch.param))
        if(length(batch.param)) {
          dbg("[correctedX] length.batch.param[1] = ", batch.param[1])
        }
        has.batchparam <- (length(batch.param)> 0 && batch.param[1]!="")
        dbg("[correctedX] has.batchparam = ", has.batchparam)
        if(m %in% c("limma","ComBat") && !has.batchparam) {
          shinyalert::shinyalert(
            type = "error",
            text = "Limma and ComBat need a valid batch parameter. Please try using SVA, RUV3 or NNM instead."
          )
        }
        
        pgx.showSmallModal("Computing batch correction. Please wait...")
        shiny::withProgress(
          message = "Computing batch correction...", value = 0.1, {
            cX <- X
            if(m %in% c("limma","ComBat") && has.batchparam) {
              mod1 <- model.matrix( ~pheno)
              batch <- samples[,batch.param]
              batch[is.na(batch)] <- "NA"  ## na not allowed
              if(m == "ComBat"  ) {
                res <- try( sva::ComBat(X, batch = batch, mod = mod1) )
                if(!"try-error" %in% class(res))  cX <- res
              }
              if(m == "limma"  ) {
                res <- try(limma::removeBatchEffect(X, batch = batch, design = mod1))
                if(!"try-error" %in% class(res))  cX <- res                
              }
            }
            if(m == "SVA")  cX <- playbase::svaCorrect(X, pheno)
            if(m == "RUV3") cX <- playbase::ruvCorrect(X, pheno)
            if(m == "NNM")  {
              dbg("[correctedX:gx.nnncorrect2] sum(is.na(X)) = ", sum(is.na(X)))
              cX <- playbase::gx.nnmcorrect2(X, pheno)$X
            }
          })        
        shiny::removeModal()
        
        cX
      }, ignoreInit = FALSE, ignoreNULL = FALSE)
      
      ## return object
      correctedCounts <- reactive({
        X <- correctedX()
        pmax( 2**X - 1, 0)
      })

      ##------------------------------------------------------------------
      ## Event listeners
      ##------------------------------------------------------------------
      
      run_outlier_methods <- eventReactive({
        list( imputedX() )
      }, {
        dbg("[outlier_server:run_outlier_methods] reacted!")       
        X <- imputedX()
        shiny::validate( shiny::need(!is.null(X), "no data. please upload."))
        shiny::validate( shiny::need(!is.null(nrow(X)), "no data. please upload."))

        X <- head(X[order(-matrixStats::rowSds(X)),],1000)
        out <- playbase::detectOutlierSamples(X, plot=FALSE, y=NULL)
        
        nb <- min(30, dim(X)/5)
        scaledX <- t(scale(t(scale(t(X), scale=FALSE))))
        corX <- cor(t(scaledX))
        
        ## standard dim reduction methods
        pos <- list()        
##        pos[['tsne']] <- Rtsne::Rtsne(scaledX, check_duplicates=FALSE, perplexity=nb)$Y

        dbg("[run_outlier_methods:irlba] dim.scaledX = ", dim(scaledX))
        pos[['pca']]  <- irlba::irlba(scaledX, nu=2, nv=0)$u
##        pos[['umap']] <- uwot::umap(scaledX, n_neighbors = ceiling(nb/2))
        for(i in 1:length(pos)) {
          rownames(pos[[i]]) <- rownames(scaledX)
          colnames(pos[[i]]) <- paste0(names(pos)[i],'_',1:2)
        }

        out$pos <- pos
        out$corX <- corX
        out
      })

      
      
      ##------------------------------------------------------------------
      ## Plot functions
      ##------------------------------------------------------------------

      plot_normalization <- function() {

##        rX <- playbase::read.as_matrix("/home/kwee/Downloads/raw_allorenteizquierdo@health.ucsd.edu_7505cdba5/raw_counts.csv")

        rX <- r_X()
        X0 <- log2( 1e-10 + rX )
        X1 <- normalizedX()

        if(input$norm_plottype == "boxplot") {
          if(ncol(X0) > 40) {
            jj <- sample(ncol(X0),40)
            ii <- rownames(X0)
            ## just downsampling for boxplots          
            if(nrow(X0) > 1000) ii <- sample(nrow(X0),1000)  
            X0 <- X0[ii,jj]
            X1 <- X1[ii,jj]
            rX <- rX[ii,jj]
          }
          
        ymax <- max( max(X0, na.rm=TRUE), max(X1, na.rm=TRUE))
          ymin <- quantile( X0[which(rX>0)], probs=0.001, na.rm = TRUE)
          
          dy <- 0.1*(ymax - ymin)
          ylim <- c( ymin - dy , ymax + dy )
          
          par(mfrow=c(1,2), mar=c(3.2,3,2,0.5), mgp=c(2.1,0.8,0) )
          boxplot( X0, main = "raw",
            ylim = ylim,
            ylab = 'expression (log2)', xlab = "samples")
          boxplot( X1, main = "normalized",
            ylim = ylim,
            ylab = '', xlab = "samples")
        }
        

      }
      

      ## missing values
      plot_missingvalues <- function() {
        X0 <- r_X()
        X1 <- imputedX()
        jj <- which(is.na(X0))
        if( isolate(input$zero_as_na) ) {
          jj <- which( is.na(X0) | X0 == 0 )
        }
        q999 <- quantile(X1, probs=0.999)[1]
        X1[X1>q999] <- NA
        h <- hist( X1, breaks=60, plot=FALSE)
        hh <- h$breaks

        par(mfrow=c(1,2), mar=c(3.2,3.2,0.8,0.5), mgp=c(2.2,0.85,0) )

        if(length(jj) > 0) {
          hist( X1[-jj], breaks=hh, main="", xlab="expression (log2CPM)" )
          hist( X1[jj], breaks=hh, add=TRUE, col='red')
        } else {
          hist( X1, breaks=hh, main="", xlab="expression (log2CPM)" )
        }

        ## NA heatmap
        par(mar=c(1.8, 1.8, 1, 2), mgp=c(2.5,0.85,0) )        
        X2 <- 1 * is.na(X0)
        if(input$zero_as_na) X2[X0==0] <- 1
        X2 <- head(X2[order(-apply(X2,1,sd)),],1000)
        playbase::gx.imagemap(X2, cex = -1)
      }

      ## sample outlier PCA plot
      plot.outlierPCA <- function(pos, z, z0, shownames) {

        is.outlier <- (z > z0)       
        col1 <- "grey70"
        ##col1 <- res.outliers$dbscan$cluster + 1
        cex1 <- cut( nrow(pos), breaks = c(0,40,100,250,1000,999999),
          c(1, 0.85, 0.7, 0.55, 0.4))
        cex1 <- 3 * as.numeric(as.character(cex1))
        pos <- playbase::uscale(pos)

        ## How about plotly??
        plot( pos, col=col1, cex = 0.8*cex1, pch = 20,
          xlim = c(-0.1,1.1), ylim = c(-0.1,1.1),
          xlab = "PC1", ylab = "PC2", main="outliers")

        if(shownames) {
          pos1 <- pos
          j <- which(is.outlier)
          if(length(j)) pos1 <- pos[-j,,drop=FALSE]
          text( pos1, rownames(pos1), cex=0.85, offset=0.8, pos = 1:4)
        }

        if(any(is.outlier)) {
          j <- which(is.outlier)
          points( pos[j,,drop=FALSE], col='red', cex = 0.8*cex1, lwd=3, pch=1 )
          outlier.name <- rownames(pos)[j]
          text( pos[j,1], pos[j,2], outlier.name, cex=1.0, offset=0.8, pos = 1:4)
        }        
      }
      
      
      ## sample outlier scores
      plot_outliers <- function() {

        res <- run_outlier_methods()

        z0 <- as.numeric(input$outlier_threshold)
        zscore <- res$z.outlier
        Z <- res$Z

        res2 <- cleanX()
        ##  pos <- res2$pos
        pos <- res$pos[['pca']]
        
        plottype <- input$outlier_plottype
        plottype = "score.pca"
        if(plottype == "score.pca") {
          par(mfrow=c(1,2), mar=c(3.2,3,2,0.5), mgp=c(2.1,0.8,0) )
          barplot( zscore, main = 'outlier score',
            ylim = c(0,max(7,1.2*max(Z))), ylab = "z-score")
          abline(h = z0, lty=3, lwd=1.5, col='red')
          plot.outlierPCA(pos, zscore, z0, input$outlier_shownames)          
        }
        
        if(plottype == "all.scores") {
          par(mfrow=c(2,3), mar=c(3,3,2,1))
          playbase::plotOutlierScores(
            res.outliers = res,
            z.threshold = z0, par=TRUE)
        }
        
      }      

      plot_correction <- function() {

        out.res <- run_outlier_methods()
        res0 <- cleanX()
        x0 <- res0$X
        x1 <- correctedX()

        dbg("[plot_correction:irlba] 1: dim.x0 = ", dim(x0))            
        dbg("[plot_correction:irlba] 1: dim.x1 = ", dim(x1))            

        ## faster
        x0 <- head(x0[order(-matrixStats::rowSds(x0,na.rm=TRUE)),],2000)
        x1 <- head(x1[order(-matrixStats::rowSds(x1,na.rm=TRUE)),],2000)

        kk <- intersect(colnames(x0),colnames(x1))
        x0 <- x0[,kk,drop = FALSE]
        x1 <- x1[,kk,drop = FALSE]

        dbg("[plot_correction:irlba] 2: dim.x0 = ", dim(x0))            
        dbg("[plot_correction:irlba] 2: dim.x1 = ", dim(x1))            
        
        plottype <- input$bec_plottype
        shownames <- input$outlier_shownames
        ## plottype <- "pca"
        if(plottype %in% c("pca","tsne","umap") ) {

          ##pos0 <- res0$pos
          if( plottype == "pca" ) {
            pos0 <- out.res$pos[['pca']][colnames(x0),]
            ## pos0 <- irlba::irlba( x0, nv=2 )$v
            dbg("[plot_correction:irlba] dim.x1 = ", dim(x1))            
            pos1 <- irlba::irlba( x1, nv=2 )$v            
          }

          if( plottype == "tsne" ) {
            p <- ceiling(min(30, dim(x0)/4))
            pos0 <- Rtsne::Rtsne( t(x0), perplexity = p )$Y
            pos1 <- Rtsne::Rtsne( t(x1), perplexity = p )$Y
          }
          
          if( plottype == "umap" ) {
            nb <- ceiling(min(15, dim(x0)/4))
            pos0 <- uwot::umap( t(x0), n_neighbors = nb )
            pos1 <- uwot::umap( t(x1), n_neighbors = nb )
          }
          
          rownames(pos0) <- colnames(x0)
          rownames(pos1) <- colnames(x1)            
          pos0 <- pos0[colnames(x0),]
          pos1 <- pos1[colnames(x0),]
          
          ##pheno <- r_contrasts()[,1]
          pheno <- playbase::contrasts2pheno( r_contrasts(), r_samples() )
          pheno <- pheno[rownames(pos0)]
          col1  <- factor(pheno)
          cex1 <- cut( nrow(pos1), breaks = c(0,40,100,250,1000,999999),
                      c(1, 0.85, 0.7, 0.55, 0.4))
          cex1 <- 2.7 * as.numeric(as.character(cex1))
          
          par( mfrow=c(1,2), mar=c(3.2,3,2,0.5), mgp=c(2.1,0.8,0) )          
          plot( pos0, col=col1, pch=20, cex = 0.9*cex1, main="before",
               xlab="PC1", ylab="PC2")
          plot( pos1, col=col1, pch=20, cex = 0.9*cex1, main="after",
               xlab="PC1", ylab="PC2")          
        }

        if(plottype == "heatmap") {
          plot.heatmap <- function(x, main) {
            x <- head(x[order(-apply(x, 1, sd)), ], 1000)
            x <- x - rowMeans(x)
            x <- abs(x)**0.5 * sign(x)
            playbase::gx.imagemap(x, main = main, cex.main = 1.4, cex = 0)
            mtext("samples", 1, line = 0.5)
            mtext("genes", 2, line = 0.5)
          }
          par( mfrow=c(1,2), mar=c(2,2,2,1), mgp=c(2.1,0.8,0) )          
          plot.heatmap( x0, main = "before")
          plot.heatmap( x1, main = "after" )
        }
        
      }
    
      ##------------------------------------------------------------------
      ## Plot modules
      ##------------------------------------------------------------------

      PlotModuleServer(
        "plot1",
        plotlib = "base",
        func = plot_normalization,
##      func2 = plot.RENDER,
##      csvFunc = plot_data,
        res = c(75, 120),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot2",
        plotlib = "base",
        func = plot_missingvalues,
        res = c(75, 120),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot3",
        plotlib = "base",
        func = plot_outliers,
##      func2 = plot.RENDER,
##      csvFunc = plot_data,
        res = c(75, 120),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot4",
        plotlib = "base",
        func = plot_correction,
        ## func2 = plot.RENDER,
        ## csvFunc = plot_data,
        res = c(75, 120),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )
      

      return( correctedCounts ) ## pointing to reactive
    } ## end-of-server
  )
}
