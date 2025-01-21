##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

DeepNetBoard <- function(id, pgx, board_observers = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 700 ## full height of page
    rowH1 <- 250 ## row 1 height
    rowH2 <- 440 ## row 2 height

    infotext <- tspan("<b>Weighted gene co-expression network analysis (WGCNA)</b> is a systems biology method for describing the correlation patterns among genes across microarray samples. Weighted correlation network analysis can be used for finding clusters (modules) of highly correlated genes, for summarizing such clusters using the module eigengene or an intramodular hub gene, for relating modules to one another and to external sample traits (using eigengene network methodology), and for calculating module membership measures. Correlation networks facilitate network based gene screening methods that can be used to identify candidate biomarkers or therapeutic targets.

<p>References:<br>
<ol>
<li>Langfelder, P. and Horvath, S., 2008. WGCNA: an R package for weighted correlation network analysis. BMC bioinformatics, 9(1), p.559.
<li>Zhang, B. and Horvath, S., 2005. A general framework for weighted gene co-expression network analysis. Statistical applications in genetics and molecular biology, 4(1).
</ol>", js = FALSE)

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Model training" = list(disable = c("show_conditions")),
      "Gradient vs. foldchange" = list(disable = NULL)
    )
    shiny::observeEvent(input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })

    
    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    infotext <-
      '<center><iframe width="1120" height="630" src="https://www.youtube.com/embed/rRIRMW_RRS4"
        title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write;
        encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>'

    shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>WGCNA Analysis Board</strong>"),
        shiny::HTML(infotext),
        size = "xl",
        easyClose = TRUE
      ))
    })

    update <- reactiveVal(0)
    
    shiny::observeEvent(
      #list(pgx$samples, pgx$multitarget)
      list(pgx$samples)
     , {
      phenotypes <- playbase::pgx.getCategoricalPhenotypes(pgx$samples)
      shiny::updateSelectInput(session, "selected_pheno", choices = phenotypes,
                               ## options = list(maxItems=2),
                               selected = phenotypes[1] )
      update( update() + 1)
    })
    
    
    shiny::observeEvent( input$selected_pheno, {
      shiny::req(input$selected_pheno)
      dbg("[DeepNetBoard] input$selected_pheno = ", input$selected_pheno)
      conditions <- sort(unique(pgx$samples[, input$selected_pheno]))
      dbg("[DeepNetBoard] conditions = ", conditions)
      shiny::updateSelectInput(session, "show_conditions", choices = conditions,
                               selected = head(conditions,3) )      
    })
    
    warned <- TRUE
    observeEvent( input$step, {
      optim <- "adam"
      #optim <- input$optim
      net()$fit( niter=1, optim=optim)
      warned <<- FALSE
      update( update() + 1)
    })
    
    observeEvent( input$step20, {
      pgx.showSmallModal(paste("Fitting model 20 steps... please wait"))
      optim <- "adam"
      #optim <- input$optim
      net()$fit( niter=20, optim=optim)
      shiny::removeModal()
      warned <<- FALSE
      update( update() + 1)
    })

    observeEvent( input$step100, {
      pgx.showSmallModal(paste("Fitting model 100 steps... please wait"))
      optim <- "adam"
      #optim <- input$optim
      net()$fit( niter=100, optim=optim)
      shiny::removeModal()
      warned <<- FALSE
      update( update() + 1)
    })

    observeEvent({
      list( 
        input$latent_dim,
        input$actfun,
        input$l1regularization,
        input$l2regularization
      )
    }, {
      if(!warned) {
        shinyalert::shinyalert("","Please reset and run your network again after adjusting network options.")
        warned <<- TRUE
      }
    })


    shiny::observeEvent({
      list( pgx$X, input$addgsets)
    },{
      shiny::req(pgx$X)
      shiny::req(!is.null(input$addgsets) && length(input$addgsets))
      
      if(!all(grepl(":",rownames(pgx$X)))) {
        shiny::updateSelectInput(session, "show_datatypes", choices = NULL,
                                 selected = NULL)
        return(NULL)
      }
      ## update datatype selectinput
      datatypes <- sort(unique(sub(":.*","",rownames(pgx$X))))
      if(input$addgsets) datatypes <- c(datatypes, "gset")
      sel.datatype <- unique(c(intersect(c("px","gx"), datatypes),datatypes))
      if(length(sel.datatype)==0) sel.datatype <- datatypes[1]
      shiny::updateSelectInput(session, "show_datatypes", choices = datatypes,
                               selected = sel.datatype[1] )

    })

    ## update network diagram if model changes and reset
    update_diagram <- reactiveVal("abc")

    
    ## ================================================================================
    ## ========================== BOARD FUNCTIONS =====================================
    ## ================================================================================
    
    ## create reactive DeepNet object
    net <- shiny::eventReactive({
      list( input$selected_pheno, input$reset )
    }, {
      
      pheno <- input$selected_pheno
      shiny::req(pheno)

      dbg("[DeepNetBoard] initialize deepnet model for",pheno)
      X <- pgx$X
      if(any(is.na(X))) {
        dbg("[DeepNetBoard] imputing missing values in X")
        X <- playbase::svdImpute2(X)
      }
      y <- pgx$samples[, pheno, drop=FALSE]
      ii <- which(rowSums(is.na(y)) == 0)
      y <- y[ii,,drop=FALSE]
      X <- X[,ii]
      sdX <- matrixStats::rowSds(X, na.rm=TRUE)
      dbg("[DeepNetBoard] 1: dim(X) = ", dim(X))
      dbg("[DeepNetBoard] 1: length(y) = ", length(y))                
      xx <- playbase::mofa.split_data(X)  ## also handles single-omics

      if(input$addgsets) {
        dbg("[DeepNetBoard] adding genesets ")
        gsetx <- pgx$gsetX[,ii]
        xx <- c( xx, list(gset=gsetx))
      }
      
      if(input$augment) {
        ntime = 10
        xx <- playbase::mofa.augment(xx, ntime, z=1)
        y  <- do.call( rbind, rep(list(y),ntime))
        dbg("[DeepNetBoard] augmented X: dim(xx[[1]]) = ", dim(xx[[1]]))
        dbg("[DeepNetBoard] augmented y: dim(y) = ", dim(y))                
      }
      
      l1 = 1
      l2 = 100
      ae.wt=1e6
      ae.wt <- 3
      if( input$model == "AE")  ae.wt <- 1e3
      if( input$model == "MLP") ae.wt <- 1e-3      
      loss_weights <- c(y=1/ae.wt, ae=ae.wt, l1=l1, l2=l2)
      dbg("[DeepNetBoard] loss_weights = ", loss_weights)                
      
      ## create the Neural network
      yy <- as.list(y)  ## one target for now...
      yy <- lapply(yy, as.factor)
      names(yy) <- pheno
      
      dbg("[DeepNetBoard] creating DeepNet model: ", input$model)
      dbg("[DeepNetBoard] layers = ", input$layers)                      
      latent_dim <- input$latent_dim

      if(input$layers == "mini") {
        num_layers = list(c(latent_dim), c(), c())
      } else if(input$layers == "deep") {
        num_layers = list(c(0.5, latent_dim), c(0.99, 0.5, 025), c(0.99, 0.5, 0.25))
      } else {
        num_layers = list(c(0.5,latent_dim), c(0.5), c(0.5))
      }
      
      net <- playbase::MultiOmicsSAE$new(
        xx,
        yy,
        model = "MT",
        num_layers = num_layers,
        ntop = 1000,
        loss_weights = loss_weights,
        #dropout = as.numeric(0.10 * input$dropout),
        add_noise = as.numeric(0.20 * input$addnoise),
        #actfun = input$actfun,
        actfun = "leaky",        
        use_glu = 2*input$useGLU,  ## GLU mode
        use_bn = TRUE
        # use_bn = input$useBN                
        # scale = input$scaleinput,
        # sd_weight = input$sdweight
      )
      dbg("[DeepNetBoard] finished creating model!")

      if( input$layers != update_diagram()) {
        update_diagram(input$layers)
      }
      warned <<- FALSE
      return(net)
    })

    
    phenoFC <- eventReactive({
      list(input$selected_pheno)
    } , {
      ## Calculate correspoding T-test statistics
      dbg("[DeepNetBoard] calculating FC..")      
      pheno <- input$selected_pheno
      y <- pgx$samples[,pheno]
      X <- pgx$X
      if(!all(grepl("[:]",rownames(X)))) {
        rownames(X) <- paste0("gx:",rownames(X))        
      }
      
      gsetX <- pgx$gsetX
      rownames(gsetX) <- paste0("gset:",rownames(gsetX))
      X <- rbind(X, gsetX)
      
      ii <- which(!is.na(y))
      y <- y[ii]
      X <- X[,ii]      
      
      dbg("[DeepNetBoard] pheno = ", pheno )
      dbg("[DeepNetBoard] y = ", y )
      
      mx <- model.matrix( ~ 0 + y)
      mx <- t(t(mx) / colSums(mx))

      dbg("[DeepNetBoard] dim.pgx$X = ", dim(X) )
      dbg("[DeepNetBoard] dim.mx = ", dim(mx) )            
      
      avgX <- X %*% mx
      F <- sapply( 1:ncol(mx), function(i) avgX[,i,drop=FALSE] - rowMeans(avgX[,-i,drop=FALSE]))
      colnames(F) <- sub("^y","",colnames(mx))
      rownames(F) <- rownames(X)
      fc <- playbase::mofa.split_data(F)
      return(fc)
    })

    
    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    plot_deepnet_diagram_server(
      "deepnet_diagram",
      net = net,  ## not reacting, only read
      pgx = pgx,  ## react on pgx
      update = update_diagram,
      watermark = WATERMARK
    )

    plot_deepnet_gradients_server(
      "deepnet_gradients",
      net = net,
      update = update,
      type = "barplot",
      conditions = reactive(input$show_conditions),
      datatype = reactive(input$show_datatypes),
      phenoFC = phenoFC,
      watermark = WATERMARK
    )

    plot_deepnet_clusters_server(
      "deepnet_clusters",
      net = net,
      update = update,
      watermark = WATERMARK
    )

    plot_deepnet_aescatter_server(
      "deepnet_aescatter",
      net = net,
      update = update,
      watermark = WATERMARK
    )

    plot_deepnet_confusionmatrix_server(
      "deepnet_confusionmatrix",
      net = net,
      update = update,
      watermark = WATERMARK
    )

    plot_deepnet_gradients_server(
      "deepnet_fcvsgrad",
      net = net,
      update = update,
      type = "scatter",
      conditions = reactive(input$show_conditions),
      datatype = reactive(input$show_datatypes),      
      phenoFC = phenoFC,
      watermark = WATERMARK
    )

    table_deepnet_gradients_server(
      "deepnet_table",
      net = net,
      pgx = pgx,
      phenoFC = phenoFC,      
      conditions = reactive(input$show_conditions),
      datatype = reactive(input$show_datatypes)
    )

    plot_deepnet_lossplot_server(
      "deepnet_lossplot",
      net = net,
      update = update,
      watermark = WATERMARK
    )

    plot_deepnet_biomarkerheatmap_server(
      "deepnet_biomarkerheatmap",
      net = net,
      update = update,
      datatypes = reactive(input$show_datatypes),
      watermark = WATERMARK
    )
    
    return(NULL)
  })
} ## end of Board
