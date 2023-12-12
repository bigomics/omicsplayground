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
        shiny::actionButton(ns("bc_compute_button"), "Batch correct", class = "run-button"),
        shiny::br(),
        shiny::br(),        
        shiny::selectInput(ns("method"), "Select method", choices="uncorrected"),
        shiny::conditionalPanel(
          "input.view == 'UMAP'",
          ns = ns,
          shiny::selectInput(ns("colorby"), "Color by", choices=NULL)
        ),
        shiny::br(), shiny::br(), shiny::br(),
        withTooltip(shiny::actionLink(ns("adv_options"),
           "Advanced options", icon = icon("cog", lib = "glyphicon")),
           "Toggle advanced options.",
           placement = "top"
           ),
        shiny::br(), shiny::br(),
        shiny::conditionalPanel(
           "input.adv_options % 2 == 1",
           ns = ns,
           shiny::checkboxInput(ns("lib_correction"), "Library correction (lib)"),
           shiny::checkboxInput(ns("cell_correction"), "Cell correction (MT/ribo/CC)"),
           shiny::checkboxInput(ns("gender_correction"), "Gender correction (M/F)"),
           shiny::hr(),
           shiny::selectInput(ns("batch_params"), "Batch covariates", choices=NULL),
           shiny::br(), shiny::br(),
           shiny::actionButton(ns("compute_button"), "Recompute", class = "run-button")
        ),
        shiny::br(), shiny::br()
      )
    ),
    shiny::mainPanel(
      width = 10,
      bslib::layout_columns(
        col_widths = 6,
        row_heights = c(2,3),
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
          ns("plot2"),
          title = "Heatmap",
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
        choices <- colnames(r_samples())
        sel <- choices[1]
        shiny::updateSelectInput(session, "colorby", choices = choices, selected=sel )

        choices2 <- sort(colnames(r_samples()))
        choices2 <- c("<none>","<all>",choices2)
        sel2 <- grep("batch", choices2, ignore.case=TRUE, value=TRUE)
        if(length(sel2)==0) sel2 <- "<none>"
        shiny::updateSelectInput(session, "batch_params", choices = choices2, selected=sel2 )
      })
      
      logX <- reactive({
        x <- playbase::counts.mergeDuplicateFeatures(r_X())
        x <- log2(x + 1e-8)
        x
      })
      
      run_all_corrections <- reactive({
        
#        samples <- pgx$samples
#        contrasts <- pgx$contrasts
#        X1 <- pgx$X
        
        samples <- r_samples()
        contrasts <- r_contrasts()
        X1 <- logX()
        
        M <- playbase::contrasts.convertToLabelMatrix( contrasts,  samples)
        vec <- playbase::createPhenoBatchVectors(
          labeled_contrasts = M,
          samples = samples,
          pmin = 0.2
        ) 

        ## Technical corrections 
        params <- c()
        if(input$lib_correction) {
          params = c(params, "lib")
        }
        if(input$gender_correction) {
          params = c(params, "gender")
        }
        if(input$cell_correction) {
          params = c(params, "mito","ribo","cellcycle")
        }
        dbg("[run_all_corrections] params = ",params)
        
        if(length(params) > 0) {
          nv <- ifelse( ncol(X1) > 40, 2, 1)          
          X1 <- playbase::removeTechnicalEffects(
            X1,
            y = vec$pheno,
            params = params,
            pmin = 0.20,
            nmin = 3,
            nv = nv
          )
        }

        ## Explicit batch corrections
        batch_params <- input$batch_params
        if(batch_params == '<none>') {
          vec$batch <- NULL
        } else if (batch_params == '<all>') {
          ## batch_params <- colnames(samples)
        } else if(batch_params %in% colnames(samples)) {
          batch_params <- intersect(batch_params, colnames(samples))
          B <- samples[,batch_params,drop=FALSE]
          B <- playbase::expandPhenoMatrix(B)
          bb <- apply(B, 1, paste0, collapse='_')
          bb <- paste0("b",as.integer(factor(bb)))
          vec$batch <- bb
        } else {
        }

        dbg("[run_all_corrections] batch_params = ",batch_params)        
        
        ## X1 <- playbase::normalizeCounts(2**X1, method='m4')
        X1 <- playbase::logCPM(pmax(2**X1 - 1,0))

        pgx.showSmallModal("Computing batch correction methods. Please wait...")
        shiny::withProgress(message = "Computing methods. Please wait...", value = 0.33, {
          
          xlist <- playbase::runBatchCorrectionMethods(X1, vec$batch, vec$pheno,
            controls=NULL, combatx = FALSE, remove.failed=TRUE) 

          incProgress( amount = 0.33)

          ## PCA is faster than UMAP
          pos <- lapply(xlist, function(x) {
            x <- head( x[order(-matrixStats::rowSds(x,na.rm=TRUE)),], 1000)
            x <- x - rowMeans(x,na.rm=TRUE)
            irlba::irlba(x)$v[,1:2]
          })

          incProgress( amount = 0.33)          
          
          res <- playbase::bc.evaluateResults(
              xlist, vec$pheno, lfc=0.2, q=0.05, pos=pos, add.sil=TRUE,
              plot=FALSE, trend=TRUE) 
        })
        shiny::removeModal()
        
        ## update selectInput
        xnames <- sort(setdiff(names(xlist), 'uncorrected'))
        choices = c("uncorrected", xnames)
        shiny::updateSelectInput(session, "method", choices = choices )

        list( xlist = xlist, pos = res$pos,
             stats = res$stats, scores = res$results )
      })
      

      output_plot1 <- function() {
        shiny::req(input$colorby, input$method)        
        res <- run_all_corrections()
        X0 <- logX()
        samples <- r_samples()
        pheno <- samples[,input$colorby]
        xlist <- res$xlist
        if(input$method != 'uncorrected') {
          kk <- c('uncorrected',input$method)        
          xlist <- xlist[kk]
        }
        
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
        res <- run_all_corrections()
        X0 <- logX()
        samples <- r_samples()
        xlist <- res$xlist
        if(input$method != 'uncorrected') {
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
          scores = res$scores,
          nmax = 1000
        ) 
      }
      
      output_plot3 <- function() {
        res <- run_all_corrections()
        X0 <- logX()
        samples <- r_samples()
        xlist <- res$xlist
        if(input$method != 'uncorrected') {
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
        res <- run_all_corrections()
        X0 <- logX()
        samples <- r_samples()
        pheno = NULL
        xlist <- res$xlist
        if(input$method != 'uncorrected') {
          kk <- c('uncorrected',input$method)        
          xlist <- xlist[kk]
        }
        playbase::bc.plotResults(
          X = X0,
          xlist = xlist,
          pos = res$pos,
          pheno = pheno,
          samples = samples,
          type = "scores",
          scores = res$scores,
          nmax = 1000
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
      
      ## ------------------------------------------------------------
      ## Reactive return object
      ## ------------------------------------------------------------
      outobj <- shiny::eventReactive({
        list(input$compute_button, input$method, r_X(), r_samples())
      }, {
        shiny::req( r_X())
        X0 <- r_X()
        dbg("[upload_module_batchcorrect_server:outobj] input$method = ", input$method)
        list(X = X0)
      },
      ignoreInit = FALSE,
      ignoreNULL = FALSE
      )

      return(outobj) ## pointing to reactive
    } ## end-of-server
  )
}
