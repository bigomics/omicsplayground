##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## =============================================================================
## ==================== BATCHCORRECT UI/SERVER =================================
## =============================================================================


upload_module_batchcorrect_ui <- function(id, height = 720) {
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
        ##        shiny::selectInput(ns("view"), "Select view", choices=c("UMAP","heatmap","PC","hist")),
        shiny::radioButtons( ns("view"), "Select view",
          choices = c("UMAP", "heatmap", "PC", "PC2", "hist"), inline=FALSE ),        

        shiny::conditionalPanel(
          "input.view == 'UMAP'",
          ns = ns,
          shiny::selectInput(ns("colorby"), "Color by", choices=NULL)
        ),
        shiny::br(),
        shiny::checkboxInput(ns("technical_correction"), "Technical correction (lib)"),
        shiny::checkboxInput(ns("cell_correction"), "Cell correction (MT/ribo/CC)"),
        shiny::checkboxInput(ns("gender_correction"), "Gender correction (M/F)")                
      )
    ),
    shiny::mainPanel(
      shiny::plotOutput(ns("canvas"), width = "100%", height = height) %>% bigLoaders::useSpinner(),
      width = 10
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

        params <- c()
        if(input$technical_correction) {
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

        ## X1 <- playbase::normalizeCounts(2**X1, method='m4')
        X1 <- playbase::logCPM(pmax(2**X1 - 1,0))
        
        shiny::withProgress(message = "Computing methods. Please wait...", value = 0.33, {
        
          xlist <- playbase::runBatchCorrectionMethods(X1, vec$batch, vec$pheno,
            controls=NULL, combatx = FALSE, remove.failed=TRUE) 

          incProgress( amount = 0.33)

          pos <- lapply(xlist, function(x) {
            x <- head( x[order(-matrixStats::rowSds(x,na.rm=TRUE)),], 1000)
            x <- x - rowMeans(x,na.rm=TRUE)
            irlba::irlba(x)$v[,1:2]
          })

          incProgress( amount = 0.33)          
          
          res <- playbase::bc.evaluateResults(
            xlist, vec$pheno, lfc=0.2, q=0.05, pos=pos, add.sil=TRUE, plot=FALSE, trend=TRUE) 
        })
        
        list( xlist = xlist, pos = res$pos, stats = res$stats )
      })
      
      
      output$canvas <- shiny::renderPlot({

        shiny::req(input$colorby)
        
        res <- run_all_corrections()
        X0 <- logX()
        samples <- r_samples()
        pheno <- samples[,input$colorby]
        viewtype <- input$view
        
        playbase::bc.plotResults(
          X = X0,
          xlist = res$xlist,
          pos = res$pos,
          pheno = pheno,
          samples = samples,
          type = viewtype,
          nmax = 1000
        ) 

      })


      ## ------------------------------------------------------------
      ## Reactive return object
      ## ------------------------------------------------------------
      outobj <- shiny::eventReactive({
        list(input$bc_compute_button, r_X(), r_samples())
      }, {
        shiny::req( r_X())
        X0 <- r_X()
        list(X = X0)
      },
      ignoreInit = FALSE,
      ignoreNULL = FALSE
      )

      return(outobj) ## pointing to reactive
    } ## end-of-server
  )
}
