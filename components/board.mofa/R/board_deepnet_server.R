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


    ## update network diagram if model changes and reset
    update <- reactiveVal(0)
    update_diagram <- reactiveVal(FALSE)
    
    ## ================================================================================
    ## ============================ OBSERVERS =========================================
    ## ================================================================================
       
    my_observers <- list()
    
    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Model training" = list(disable = c("show_conditions")),
      "Gradient vs. foldchange" = list(disable = NULL),
      "Biomarker heatmap" = list(disable = "show_conditions")
    )

    my_observers[[1]] <- shiny::observeEvent(input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })

    infotext2 <-
      '<center><iframe width="1120" height="630" src="https://www.youtube.com/embed/rRIRMW_RRS4"
        title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write;
        encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>'

    my_observers[[2]] <- shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>WGCNA Analysis Board</strong>"),
        shiny::HTML(infotext2),
        size = "xl",
        easyClose = TRUE
      ))
    })

    
    my_observers[[3]] <-  shiny::observeEvent({
      list(pgx$samples)
    } , {
      phenotypes <- playbase::pgx.getCategoricalPhenotypes(pgx$samples)
      shiny::updateSelectInput(session, "selected_pheno", choices = phenotypes,
                               ## options = list(maxItems=2),
                               selected = phenotypes[1] )
      update( update() + 1)
    })
    
    
    my_observers[[4]] <- shiny::observeEvent({
      list(input$selected_pheno, pgx$samples)
    }, {
      shiny::req(input$selected_pheno)
      shiny::req(pgx$samples)
      if(!input$selected_pheno %in% colnames(pgx$samples)) return(NULL)
      conditions <- sort(unique(pgx$samples[, input$selected_pheno]))
      shiny::updateSelectInput(session, "show_conditions", choices = conditions,
                               selected = head(conditions,3) )      
    })
    
    warned <- TRUE
    my_observers[[5]] <- observeEvent( input$step, {
      optim <- "adam"
      #optim <- input$optim
      net()$fit( niter=1, optim=optim)
      warned <<- FALSE
      update( update() + 1)
    })
    
    my_observers[[6]] <- observeEvent( input$step20, {
      pgx.showSmallModal(paste("Fitting model 20 steps... please wait"))
      optim <- "adam"
      #optim <- input$optim
      net()$fit( niter=20, optim=optim)
      shiny::removeModal()
      warned <<- FALSE
      update( update() + 1)
    })
    
    my_observers[[7]] <- observeEvent( input$step100, {
      pgx.showSmallModal(paste("Fitting model 100 steps... please wait"))
      optim <- "adam"
      #optim <- input$optim
      net()$fit( niter=100, optim=optim)
      shiny::removeModal()
      warned <<- FALSE
      update( update() + 1)
    })

    my_observers[[8]] <- observeEvent({
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

    my_observers[[9]] <- shiny::observeEvent({
      list(pgx$X, input$addgsets)
    },{
      
      shiny::req(pgx$X)
      shiny::req(!is.null(input$addgsets) && length(input$addgsets))
      
      if(!all(grepl(":",rownames(pgx$X)))) {
        dbg("[DeepNetBoard] 0: SINGLE OMICS???")
        shiny::updateSelectInput(
          session, "show_datatypes", choices = "gx", selected = "gx")
        return(NULL)
      }

      ## update datatype selectinput
      sel.datatype <- input$show_datatypes

      datatypes <- sort(unique(sub(":.*","",rownames(pgx$X))))
      if(input$addgsets) {
        datatypes <- unique(c(datatypes, "GSET"))
        sel.datatype <- unique(c(sel.datatype, "GSET"))
      }
      sel.datatype <- intersect(sel.datatype, datatypes)

      if(length(sel.datatype)==0) sel.datatype <- datatypes
      shiny::updateSelectInput(session, "show_datatypes", choices = datatypes,
        selected = sel.datatype )

    })

    my_observers[[10]] <- shiny::observeEvent({
      ##list(input$addgsets, input$show_datatypes, input$layers)
      list( input$selected_pheno, input$reset )
    },{
      update_diagram(TRUE)
    })
    
    ## add to list global of observers. suspend by default.
    my_observers <- my_observers[!sapply(my_observers,is.null)]
    # lapply( my_observers, function(b) b$suspend() )
    if(!is.null(board_observers)) board_observers[[id]] <- my_observers
    
    ## ================================================================================
    ## ========================== BOARD FUNCTIONS =====================================
    ## ================================================================================
    
    ## create reactive DeepNet object
    net <- shiny::eventReactive({
      list( input$selected_pheno, input$reset )
    }, {
      shiny::req(input$selected_pheno)
      shiny::req(input$show_datatypes)      

      pheno <- input$selected_pheno
      X <- pgx$X
      if(any(is.na(X))) {
        info("[DeepNetBoard] imputing missing values in X")
        X <- playbase::svdImpute2(X)
      }
      y <- pgx$samples[, pheno, drop=FALSE]
      ii <- which(rowSums(is.na(y)) == 0)
      y <- y[ii,,drop=FALSE]
      X <- X[,ii]
      sdX <- matrixStats::rowSds(X, na.rm=TRUE)
      xx <- playbase::mofa.split_data(X)  ## also handles single-omics
      
      ## Subset selection of datatypes for deepnet
      dt <- isolate(input$show_datatypes)
      if(length(dt)) {
        xx <- xx[setdiff(dt,"GSET")]
      }
      
      if(input$addgsets) {
        info("[DeepNetBoard] adding genesets ")
        gsetx <- pgx$gsetX[,ii]
        xx <- c( xx, list(GSET=gsetx))
        update_diagram(TRUE)        
      }
      
      if(input$augment) {
        ntime = 10
        xx <- playbase::mofa.augment(xx, ntime, z=1)
        y  <- do.call( rbind, rep(list(y),ntime))
      }
      
      l1 = 1
      l2 = 100
      ae.wt=3
      if( input$model == "AE")  ae.wt <- 1e3
      if( input$model == "MLP") ae.wt <- 1e-3      
      loss_weights <- c(y=1/ae.wt, ae=ae.wt, l1=l1, l2=l2)
      
      ## create the Neural network
      yy <- as.list(y)  ## one target for now...
      yy <- lapply(yy, as.factor)
      names(yy) <- pheno
      
      info("[DeepNetBoard] creating DeepNet model: ", input$model)
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
        add_noise = as.numeric(1.2 * input$addnoise),
        #actfun = input$actfun,
        actfun = "leaky",        
        use_glu = 2*input$useGLU,  ## GLU mode
        use_bn = TRUE
        # use_bn = input$useBN                
        # scale = input$scaleinput,
        # sd_weight = input$sdweight
      )

      info("[DeepNetBoard] finished creating model!")
      ## update_diagram(TRUE)      ## need update?
      warned <<- FALSE
      return(net)
    })

    
    phenoFC <- eventReactive({
      list(input$selected_pheno)
    } , {
      ## Calculate correspoding T-test statistics
      pheno <- input$selected_pheno
      y <- pgx$samples[,pheno]
      X <- pgx$X
      if(!all(grepl("[:]",rownames(X)))) {
        rownames(X) <- paste0("gx:",rownames(X))        
      }      
      gsetX <- pgx$gsetX
      rownames(gsetX) <- paste0("GSET:",rownames(gsetX))
      X <- rbind(X, gsetX)
      
      ii <- which(!is.na(y))
      y <- y[ii]
      X <- X[,ii]      
      
      mx <- model.matrix( ~ 0 + y)
      mx <- t(t(mx) / colSums(mx))
      
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
      pgx = pgx,
      ntop = c(20,32),
      rmar = c(0,40),
      show_legend = c(0,1),
      add_annot = c(0,1),
      update = update,
      watermark = WATERMARK
    )

    plot_deepnet_biomarkerheatmap_server(
      "deepnet_bigheatmap",
      net = net,
      pgx = pgx,
      ntop = c(50,50),      
      plot.res = c(110, 110),
      add_annot = c(1,1),
      show_legend = c(1,1),
      rmar = c(40,40),
      update = update,
      datatypes = reactive(input$show_datatypes),      
      watermark = WATERMARK
    )
    
    return(NULL)
  })
} ## end of Board
