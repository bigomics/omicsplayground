##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## =============================================================================
## ==================== BATCHCORRECT UI/SERVER =================================
## =============================================================================


upload_module_batchcorrect_ui <- function(id, height = "100%") {
  ns <- shiny::NS(id)

  bc_info <- "Batch correction can clean your data from 'unwanted variables'. Please specify your parameters of interest.\n"

  shiny::sidebarLayout(
    shiny::sidebarPanel(
      width = 2,
      shiny::tagList(
        shiny::p(bc_info),
        shiny::br(),
        shiny::actionButton(ns("skip_button"), "Skip",
          class = "btn-sm btn-outline-success me-1", width='100%'),
        br(), br(),
        shiny::selectizeInput(ns("method"), NULL,
          choices = c("Select method" = "") ),
        shiny::actionButton(ns("use_button"), "Use batch correction",
          class = "btn-sm btn-outline-primary me-1", width='100%'),        
        shiny::br(), shiny::br(), shiny::br(), shiny::br(),
        withTooltip(shiny::actionLink(ns("adv_options"),
           "Advanced options", icon = icon("cog", lib = "glyphicon")),
           "Toggle advanced options.",
           placement = "top"
           ),
        shiny::br(), shiny::br(),
        shiny::conditionalPanel(
           "input.adv_options % 2 == 1",
           ns = ns,
           shiny::selectizeInput(ns("tech_params"), "Technical covariates:",
             choices = c("<lib>","<mito>","<ribo>","<cell>","<gender>"),
             multiple=TRUE),           
           shiny::selectizeInput(ns("batch_params"), "Batch covariates:",
             choices = NULL, multiple=TRUE),
           shiny::actionButton(ns("recompute_button"), "Recompute",
             class="btn-sm", width='100%'),
        ),
        shiny::br(), shiny::br()
      )
    ),
    shiny::mainPanel(
      width = 10,
      bslib::layout_columns(
        col_widths = 6,
        row_heights = c(3,3),
        height = "calc(100vh - 200px)",
        heights_equal = "row",
        ##  shiny::plotOutput(ns("canvas"), width = "100%", height = height) %>% bigLoaders::useSpinner(),
        PlotModuleUI(
          ns("plot1"),
          title = "Clustering",
#          info.text = info.text,
#          caption = caption,
          options = NULL
        ),
        PlotModuleUI(
          ns("plot3"),
          title = "PC components",
#          info.text = info.text,
#          caption = caption,
          options = NULL
        ),
        PlotModuleUI(
          ns("plot2"),
          title = "Heatmap",
#          info.text = info.text,
#          caption = caption,
          options = NULL
        ),
        PlotModuleUI(
          ns("plot4"),
          title = "Statistics and score",
#          info.text = info.text,
#          caption = caption,
          options = NULL
        )        
      )
    )
  )
}


upload_module_batchcorrect_server <- function(id, r_X, r_samples, r_contrasts,
                                              is.count = FALSE, height = 720) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      
      observeEvent( r_samples(), {
        ## choices1 = c("<lib>", "<mito>", "<ribo>", "<cellcycle>", "<gender>","<none>")
        ## shiny::updateSelectInput(session, "tech_params", choices = choices1 )
        ## choices2 = c(colnames(r_samples()))
        ## shiny::updateSelectInput(session, "batch_params", choices = choices2)        
      })

      logX <- eventReactive({
        list( input$recompute_button, r_X(), r_samples(), r_contrasts() )
      } , {
        samples <- r_samples()
        contrasts <- r_contrasts()
        dbg("[batchcorrect_server:logX] reacted!")        
        x <- playbase::counts.mergeDuplicateFeatures(r_X())
        x <- limma::normalizeQuantiles(playbase::logCPM(x))
        x
      })

      analyze_batch_effects <- reactive({

        dbg("[batchcorrect_server:analyze_batch_effects] 0 :")
        dbg("[batchcorrect_server:analyze_batch_effects] 0 : dim.X = ", dim(r_X()))
        dbg("[batchcorrect_server:analyze_batch_effects] 0 : dim.contrasts = ", dim(r_contrasts()))
        dbg("[batchcorrect_server:analyze_batch_effects] 0 : dim.samples = ", dim(r_samples()))
        
        ##shiny::req(r_contrasts(), r_samples(), r_X())
        shiny::req( dim(r_contrasts()), dim(r_samples()), dim(r_X()))

        dbg("[batchcorrect_server:analyze_batch_effects] 1 :")        
        
        samples <- r_samples()
        contrasts <- r_contrasts()
        ##X <- r_X()
        X <- logX()
        
        dbg("[batchcorrect_server:analyze_batch_effects] 2 :")
        
        ## M <- contrasts
        ## if( input$model_params == "<auto>") {
        ##   M <- playbase::contrasts.convertToLabelMatrix( contrasts,  samples)
        ##   pheno <- playbase::contrasts2pheno(M, samples) 
        ## } else {
        ##   kk <- intersect(input$model_params, colnames(samples))
        ##   M <- samples[,kk,drop=FALSE]
        ##   pheno <- apply( 1*playbase::expandPhenoMatrix(M), 1, paste, collapse='')
        ## }

        dbg("[batchcorrect_server:analyze_batch_effects] 3 :")
        
        bc <- playbase::detectBatchEffects(
          X = X,
          samples = samples,
          pheno = NULL,
          contrasts = contrasts,
          ##params = c("statistical","technical","pca"),
          params = c("statistical","technical"),          
          p.pca = 0.5,
          p.pheno = 0.05,
          k.pca = 10,
          nv = 2,
          xrank = NULL
        )                

        dbg("[batchcorrect_server:analyze_batch_effects] 4 :")
        
        choices2 = c("<none>",sort(bc$params$statistical))
        shiny::updateSelectInput(session, "batch_params", choices = choices2, sel = choices2[-1])

        tparams <- sort(unique(sub("[.].*","",bc$params$technical)))
        choices3 <- c("<none>",paste0("<",tparams,">"))
        shiny::updateSelectInput(session, "tech_params", choices = choices3, sel = choices3[-1])

        dbg("[batchcorrect_server:analyze_batch_effects] 5 :")
        
        bc$choices <- list(
          batch_params = choices2[-1],
          tech_params  = choices3[-1]          
        )
        
        bc
      })
      
      run_methods <- eventReactive({
        list( logX(), input$recompute_button )
      } , {
        
        dbg("[batchcorrect_server:run_methods] reacted!")
        
        samples <- r_samples()
        contrasts <- r_contrasts()
        X1 <- logX()

        dbg("[batchcorrect_server:run_methods] 1 :")
        
        pgx.showSmallModal("Computing batch correction methods. Please wait...")

        shiny::withProgress(message = "Computing methods. Please wait...", value = 0.1, {

          shiny::incProgress( amount = 0.2, "Analyzing for batch effects..." )          
          bc <- analyze_batch_effects()        
          shiny::req(bc)
          
          tparams <- input$tech_params
          bparams <- input$batch_params

          dbg("[batchcorrect_server:run_methods] isnull.bparams = ", is.null(bparams))
          dbg("[batchcorrect_server:run_methods] isnull.tparams = ", is.null(tparams))
          dbg("[batchcorrect_server:run_methods] bparams=='' = ", bparams=='')
          dbg("[batchcorrect_server:run_methods] tparams=='' = ", tparams=='')
          dbg("[batchcorrect_server:run_methods] len.bparams = ", length(bparams))
          dbg("[batchcorrect_server:run_methods] len.tparams = ", length(tparams))
          if(is.null(tparams) || length(tparams)==0) {
            tparams <- bc$choices$tech_params
          }
          if(is.null(bparams) || length(bparams)==0) {
            bparams <- bc$choices$batch_params
          }
          dbg("[batchcorrect_server:run_methods] 2: bparams = ", bparams)
          dbg("[batchcorrect_server:run_methods] 2: tparams = ", tparams)

          batch <- bc$batch          
          M <- bc$covariates
          mparams <- sub("=.*","",colnames(M))
          sel <- c()
          if(length(tparams)==0 || '<none>' %in% tparams) {
            ##
          } else {
            if("<lib>" %in% tparams) sel <- c(sel, grep("^lib[.]",mparams))
            if("<ribo>" %in% tparams) sel <- c(sel, grep("^ribo[.]",mparams))
            if("<mito>" %in% tparams) sel <- c(sel, grep("^mito[.]",mparams))
            if("<gender>" %in% tparams) sel <- c(sel, grep("^gender[.]",mparams))
            if("<cellcycle>" %in% tparams) sel <- c(sel, grep("^cellcycle[.]",mparams))
            if("<pca>" %in% tparams) sel <- c(sel, grep("^pca[.]",mparams))
          }
          
          if(length(bparams)==0 || '<none>' %in% bparams) {
            ##
          } else {
            sel1 <- which(mparams %in% bparams)
            sel <- c(sel, sel1)
          }

          ## create batch vector
          batch <- NULL
          if(length(sel)>0) {
            M1 <- M[,sel, drop=FALSE]
            batch <- playbase::samples2pheno(M1)
          }
          
          shiny::incProgress( amount = 0.2, "Running correction methods...")

          methods <- c("uncorrected","ComBat", "limma","superBC",
            "PCA","RUV","SVA","NNM")
          methods <- c("uncorrected","ComBat", "limma","RUV","SVA","NNM")
          xlist <- playbase::runBatchCorrectionMethods(
            X = X1,
            batch = batch,
            y = bc$pheno,
            controls = NULL,
            methods = methods,
            combatx = FALSE,
            ntop = 2000,
            sc = FALSE,
            remove.failed=TRUE)         

          incProgress( amount = 0.2, "Computing clustering...")

          ## PCA is faster than UMAP
          pos <- lapply(xlist, function(x) {
            x <- head( x[order(-matrixStats::rowSds(x,na.rm=TRUE)),], 1000)
            x <- x - rowMeans(x,na.rm=TRUE)
            irlba::irlba(x)$v[,1:2]
          })

          incProgress( amount = 0.2, "Evaluating results...")          
          
          res <- playbase::bc.evaluateResults(
            xlist,
            pheno = bc$pheno,
            lfc = 0.2,
            q = 0.05,
            pos = pos,
            add.sil = TRUE,
            plot = FALSE,
            trend = TRUE
          )
          
        })
        shiny::removeModal()
        
        ## update selectInput with methods
        xnames <- sort(setdiff(names(xlist), 'uncorrected'))
        choices = c("uncorrected", xnames)
        choices = c("Select method" = "", xnames)
        shiny::updateSelectInput(session, "method", choices = choices )

        list(
          xlist = xlist,
          pos = res$pos,
          scores = res$scores,
          pheno = bc$pheno
        )
      })


      ##      correctedX <- shiny::eventReactive({
      ##        input$use_button
      ##      }, {
      correctedX <- reactive({      
        dbg("[upload_module_batchcorrect_server:correctedX] reacted!")
        X <- logX()
        this.method <- input$method        
        dbg("[upload_module_batchcorrect_server:correctedX] this.method = ", this.method)
        bc <- analyze_batch_effects()        
        xlist <- playbase::runBatchCorrectionMethods(
          X = X,
          batch = bc$batch,
          y = bc$pheno,
          controls = NULL,
          methods = this.method,
          combatx = FALSE,
          ntop = Inf,
          sc = FALSE,
          remove.failed=TRUE)         
        xlist[[this.method]]
      })

      ## ------------------------------------------------------------
      ## Reactive return object
      ## ------------------------------------------------------------
      outobj <- shiny::reactive({
        shiny::req(correctedX())
        X0 <- correctedX()
        list(X = X0)
      })
            
      output_plot1 <- function() {
        res <- run_methods()
        shiny::req(res)
        X0 <- isolate(logX())
        samples <- isolate(r_samples())
        pheno <- res$pheno
        xlist <- res$xlist
        if(input$method %in% names(xlist) && input$method != 'uncorrected') {
          kk <- c('uncorrected',input$method)        
          xlist <- xlist[kk]
        }
        dbg("[plot1] table.pheno = ", table(pheno))
        
        playbase::bc.plotResults(
          X = X0,
          xlist = xlist,
          pos = res$pos,
          pheno = pheno,
          samples = samples,
          type = "umap",
          scores = NULL,
          nmax = 1000
          )
      }
      
      output_plot2 <- function() {
        res <- run_methods()
        shiny::req(res)
        X0 <- isolate(logX())
        samples <- isolate(r_samples())
        xlist <- res$xlist
        if(input$method %in% names(xlist) &&  input$method != 'uncorrected') {
          kk <- c('uncorrected',input$method)        
          xlist <- xlist[kk]
        }
        
        playbase::bc.plotResults(
          X = X0,
          xlist = xlist,
          pos = res$pos,
          pheno = NULL,
          samples = samples,
          type = "heatmap",
          scores = NULL,
          nmax = 1000
        ) 
      }
      
      output_plot3 <- function() {
        res <- run_methods()
        X0 <- isolate(logX())
        samples <- isolate(r_samples())
        xlist <- res$xlist
        if(input$method %in% names(xlist) &&  input$method != 'uncorrected') {
          kk <- c('uncorrected',input$method)        
          xlist <- xlist[kk]
        }

        playbase::bc.plotResults(
          X = X0,
          xlist = xlist,
          pos = res$pos,
          pheno = NULL,
          samples = samples,
          type = "pc",
          scores = NULL,
          nmax = 1000
        ) 
      }

      output_plot4 <- function() {
        res <- run_methods()
        shiny::req(res)
        X0 <- isolate(logX())
        samples <- isolate(r_samples())
        pheno = NULL
        xlist <- res$xlist
        sel <- grep("^r[.]g",colnames(res$scores),invert=TRUE)
        scores <- res$scores[,sel]
        pos = res$pos
        
        dbg("[output_plot4] rownames.scores = ",rownames(scores))
        dbg("[output_plot4] input$method = ",input$method)
        
        if(input$method %in% rownames(scores) && input$method != 'uncorrected') {
          kk <- c('uncorrected',input$method)        
          scores <- scores[kk,]
          xlist <- xlist[kk]
          pos <- pos[kk]          
        }
        
        playbase::bc.plotResults(
          X = X0,
          xlist = xlist,
          pos = pos,
          pheno = pheno,
          samples = samples,
          type = "scores",
          scores = scores,
          nmax = 1000,
          ncol = 3
        ) 
      }


      PlotModuleServer(
        "plot1",
        plotlib = "base",
        func = output_plot1,
##        func2 = plot.RENDER,
##        csvFunc = plot_data,
        res = c(70, 110),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot2",
        plotlib = "base",
        func = output_plot2,
##        func2 = plot.RENDER,
##        csvFunc = plot_data,
        res = c(70, 110),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot3",
        plotlib = "base",
        func = output_plot3,
##        func2 = plot.RENDER,
##        csvFunc = plot_data,
        res = c(70, 110),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot4",
        plotlib = "base",
        func = output_plot4,
##        func2 = plot.RENDER,
##        csvFunc = plot_data,
        res = c(70, 110),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )
      

      return(outobj) ## pointing to reactive
    } ## end-of-server
  )
}
