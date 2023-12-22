##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## =============================================================================
## ==================== BATCHCORRECT UI/SERVER =================================
## =============================================================================


upload_module_batchcorrect_ui <- function(id, height = "100%") {
  ns <- shiny::NS(id)

  bc_info <- HTML("<h4>Batch correction</h4>Batch correction can clean your data from 'unwanted variables'.\n")

  clust.options <- tagList(
    shiny::selectInput( ns('plottype'), "Type;", c("tsne","pca","heatmap"))
  )
  clust.options <- tagList(
    shiny::radioButtons( ns('plottype'), "Type;", c("tsne","pca","heatmap"))
  )

  clust.infotext =
  "Clustering of samples before ('uncorrected') and after the different batch-effect correction methods. After batch-effect correction clustering should improve. The silhouette score gives an indication of the clustering performance of the method."

  va.infotext =
  "Analysis of variables by plotting their significance in correlation with the phenotype against their significance in correlation with a principal component (PC) vector. Strong model variables are situate 'top right'. Batch effect variables with high PC correlation but low phenotype correlation are on the 'top left'. A well-designed experiment shows strong model variables in PC1, else it may be a sign of significant batch-effects."
  
  bslib::layout_columns(
    col_widths = c(2,10),
    height = "calc(100vh - 200px)",
    heights_equal = "row",
    bslib::card(
      style = 'background-color: #F7FAFD88;',
      bc_info,
      shiny::br(),
      shiny::br(),        
      shiny::selectizeInput(ns("method"), "Select method:",
        choices = c("<uncorrected>") ),
      shiny::conditionalPanel(
        "input.method != '<uncorrected>'",
        ns = ns,
        shiny::actionButton(ns("use_button"),
          label = "Apply batch correction ",
          class = "btn-sm btn-primary me-1", width='100%')
      ),
      shiny::br(),
      shiny::br(),
      withTooltip(shiny::actionLink(ns("adv_options"),
        "Advanced options", icon = icon("cog", lib = "glyphicon")),
        "Toggle advanced options.",
        placement = "top"
      ),
      shiny::conditionalPanel(
        "input.adv_options % 2 == 1",
        ns = ns,
        shiny::selectizeInput(ns("tech_params"), "Technical covariates:",
          choices = "<none>",  multiple = TRUE),           
        shiny::selectizeInput(ns("batch_params"), "Batch covariates:",
          choices = "<none>", multiple=TRUE),
        shiny::actionButton(ns("recompute_button"), "Recompute",
          class="btn-sm btn-primary mt-3", width='100%'),
        ),
      shiny::br(), shiny::br()
    ),
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
          title = "Clustering",
          info.text = clust.infotext,
          caption = clust.infotext,
          options = clust.options,
          height = c("100%","70vh")          
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
          title = "Variable analysis",
          info.text = va.infotext,
          caption = va.infotext,
          options = NULL,
          height = c("100%","70vh")
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
      
      ## also return object
      correctedX <- reactiveVal(NULL)
      observeEvent(r_X(), correctedX(r_X()))

      logX <- eventReactive({
        list( r_X() )
      } , {
        shiny::validate( shiny::need(!is.null(r_X()), "no data. please upload."))
        shiny::validate( shiny::need(!is.null(nrow(r_X())), "no data. please upload."))
        X <- playbase::counts.mergeDuplicateFeatures(r_X())
        X <- limma::normalizeQuantiles(playbase::logCPM(X))  ## standard normalization
        X
      })

      analyze_batch_effects <- reactive({

        ##shiny::req(r_contrasts(), r_samples(), r_X())
        shiny::req( dim(r_contrasts()), dim(r_samples()), dim(logX()))       
        samples <- r_samples()
        contrasts <- r_contrasts()
        X <- logX()

        shiny::validate( shiny::need(!is.null(r_X()) && nrow(r_X()),
          "No counts data. Please upload."))
        shiny::validate( shiny::need(!is.null(r_contrasts()) && nrow(r_contrasts()),
          "No contrasts defined. Please create contrasts."))
        
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

        bparams <- c()
        tparams <- c()        

        if("technical" %in% names(bc$params) && length(bc$params$technical)) {
          tparams <- sort(unique(sub("[.].*","",bc$params$technical)))
          choices3 <- c("<none>",paste0("<",tparams,">"))
          shiny::updateSelectInput(session, "tech_params", choices = choices3,
            sel = choices3[-1])
        } else {
          shiny::updateSelectInput(session, "tech_params", choices = "<none>",
            sel = "<none>")
        }
        
        if("statistical" %in% names(bc$params) && length(bc$params$statistical)) {
          bparams <- sort(bc$params$statistical)
          choices2 = c("<none>",bparams)
          shiny::updateSelectInput(session, "batch_params", choices = choices2,
            sel = choices2[-1])
        } else {
          shiny::updateSelectInput(session, "batch_params", choices = "<none>",
            sel = "<none>")
        }

        if( length(tparams) == 0 && length(bparams) == 0 ) {
          shinyalert::shinyalert(
            title = "",                        
            text = "Your data has no significant batch effects. We recommend proceeding without correction."
          )
        }
        if( length(tparams) > 0 || length(bparams) > 0 ) {
          shinyalert::shinyalert(
            type = "warning",
            title = "",                                                
            text = paste("Your data has batch effects. We recommend to use one of the correction methods.")
          )
        }
        
        bc$choices <- list(
          batch_params = bparams,
          tech_params  = tparams
        )
        
        bc
      })
      
      run_methods <- eventReactive({
        list( logX(), r_contrasts(), input$recompute_button )
      } , {
        
        dbg("[batchcorrect_server:run_methods] reacted!")
        samples <- r_samples()
        contrasts <- r_contrasts()
        X1 <- logX()

        dbg("[batchcorrect_server::run_methods] 0: dim(X1) = ",dim(X1))
                
        shiny::validate( shiny::need(!is.null(X1) && nrow(X1),
          "No counts data. Please upload."))
        shiny::validate( shiny::need(!is.null(contrasts) && ncol(contrasts),
          "No contrasts. Please create contrasts."))

        dbg("[run_methods] dim.X = ",dim(X1))
        dbg("[run_methods] dim.contrasts = ",dim(contrasts))
        
        ##pgx.showSmallModal("Computing batch correction methods. Please wait...")

        shiny::withProgress(message = "Computing methods. Please wait...", value = 0.1, {

          shiny::incProgress( amount = 0.1, "Analyzing for batch effects..." )          
          bc <- analyze_batch_effects()        
          shiny::req(bc)
          
          tparams <- input$tech_params
          bparams <- input$batch_params

          if(is.null(tparams) || length(tparams)==0) {
            tparams <- bc$choices$tech_params
          }
          if(is.null(bparams) || length(bparams)==0) {
            bparams <- bc$choices$batch_params
          }

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
          
          shiny::incProgress( amount = 0.1, "Running correction methods...")
          
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

          incProgress( amount = 0.0, "Computing clustering...")
          
          ## PCA is faster than UMAP
          pos <- list()
          t2 <- function(x) t(scale(t(scale(t(x),scale=FALSE))))
          nb <- max(2,round(min(30, dim(X1)/4)))

          incProgress( amount = 0.1, "Computing PCA clustering...")          
          pos[['pca']] <- lapply(xlist, function(x) {
            irlba::irlba(t2(x), nu=2, nv=2)$u[,1:2]
          })
          incProgress( amount = 0.1, "Computing t-SNE clustering...")          
          pos[['tsne']] <- lapply(xlist, function(x) {
            Rtsne::Rtsne(t2(x), perplexity=nb, check_duplicates=FALSE)$Y
          })
          ## incProgress( amount = 0.1, "Computing UMAP clustering...")          
          ## pos[['umap']] <- lapply(xlist, function(x) {
          ##   as.matrix(uwot::umap(t2(x), n_neighbors=nb/2))
          ## })


          
          incProgress( amount = 0.1, "Evaluating results...")          

          dbg("[run_methods] calling bc.evaluateResults")
          
          res <- playbase::bc.evaluateResults(
            xlist,
            pheno = bc$pheno,
            lfc = 0.2,
            q = 0.05,
            pos = pos[['tsne']],
            add.sil = TRUE,
            plot = FALSE,
            trend = TRUE
          )
          
        })
        ##shiny::removeModal()

        ## update selectInput with methods
        xnames <- sort(setdiff(names(xlist), 'uncorrected'))
        choices = c("<uncorrected>", xnames)
        ##choices = c("Select method" = "", xnames)
        shiny::updateSelectInput(session, "method", choices = choices )

        list(
          xlist = xlist,
          pos = pos,
          scores = res$scores,
          pheno = bc$pheno
        )
      })
      
      shiny::observeEvent({
        input$use_button
      },{
        dbg("[batchcorrect_server:correctedX] reacted!")

        this.method <- input$method

        shinyalert::shinyalert(
          text = paste("Correcting your data using",this.method,"...")
        )
        
        shiny::withProgress(message = "Correcting data...", value = 0.33, {
            bc <- analyze_batch_effects()        
            xlist <- playbase::runBatchCorrectionMethods(
              X = logX(),
              batch = bc$batch,
              y = bc$pheno,
              controls = NULL,
              methods = this.method,
              combatx = FALSE,
              ntop = Inf,
              sc = FALSE,
              remove.failed = TRUE
            )         
        })
        
        ## copy to reactive value
        corrX <- xlist[[this.method]]
        corr_counts <- pmax( 2**corrX - 1, 0)
        correctedX( corr_counts )
        ##shiny::removeModal()
        
        shinyalert::shinyalert(
          text = paste("Your data has been batch-corrected using",this.method,
                       ". You can now continue to the Compute tab."),
          immediate = TRUE
        )
        
      })
            
      output_plot1 <- function() {
        res <- run_methods()
        shiny::req(res)
        X0 <- isolate(logX())
        samples <- isolate(r_samples())
        pheno <- res$pheno
        xlist <- res$xlist
        if(input$method %in% names(xlist) &&
             !grepl("uncorrected", input$method)) {
          kk <- c('uncorrected',input$method)        
          xlist <- xlist[kk]
        }

        type = "umap"
        if(input$plottype == "heatmap") {
          type = "heatmap"
          pos <- NULL
        } else {
          type = "umap"
          layout <- input$plottype
          pos <- res$pos[[layout]]
        }
        
        playbase::bc.plotResults(
          X = X0,
          xlist = xlist,
          pos = pos,
          pheno = pheno,
          samples = samples,
          type = type,
          scores = NULL,
          nmax = 1000,
          cex = 0.65
          )
      }
      
      output_plot2 <- function() {
        bc <- analyze_batch_effects()       
        shiny::req(bc)
        par(mfrow=c(2,2), mar=c(3.3,4,2,2), mgp=c(2.2,1,0))
        playbase::bc.AnalysisPlotPCA(bc, k=1:4, par=FALSE, col='#273871')
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
          pos = NULL,
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
        sel <- grep("^r[.]g|avg.sd",colnames(res$scores),invert=TRUE)
        scores <- res$scores[,sel]
        
        if(input$method %in% rownames(scores) && input$method != 'uncorrected') {
          kk <- c('uncorrected',input$method)        
          scores <- scores[kk,]
          xlist <- xlist[kk]
        }
        
        playbase::bc.plotResults(
          X = X0,
          xlist = xlist,
          pos = NULL,
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
        res = c(65, 110),
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
        res = c(65, 110),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )
      

      return( correctedX ) ## pointing to reactive
    } ## end-of-server
  )
}
