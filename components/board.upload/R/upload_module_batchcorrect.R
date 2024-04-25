##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## =============================================================================
## ==================== BATCHCORRECT UI/SERVER =================================
## =============================================================================


upload_module_batchcorrect_ui <- function(id, height = "100%") {
  ns <- shiny::NS(id)

  bc_info <- HTML("<h4>Batch-effects analysis</h4>Batch correction can clean your data from 'unwanted variables'.\n")

  clust.infotext =
  "Clustering of samples before ('uncorrected') and after the different batch-effect correction methods. After batch-effect correction clustering should improve. The silhouette score gives an indication of the clustering performance of the method."

  covariate.info =
  "Analysis of variables by plotting their significance in correlation with the phenotype against their significance in correlation with a principal component (PC) vector. Strong model variables are situate 'top right'. Batch effect variables with high PC correlation but low phenotype correlation are on the 'top left'. A well-designed experiment shows strong model variables in PC1, else it may be a sign of significant batch-effects."

  pcc.info = "PC analysis by covariate (class). The heights of the bars correspond to the relative contribution of that covariate to the three PC, as measured by an F-test. "

  clust.options <- tagList(
    shiny::radioButtons( ns('clust.plottype'), "Type;", c("tsne","pca","heatmap"))
  )

  covariate.options <- tagList(
      shiny::radioButtons( ns('covariate.plottype'), "Type;",
                          c("Covariate plot","Correlation heatmap"))
  )

  
  bslib::layout_columns(
    col_widths = c(2,10),
    height = "calc(100vh - 200px)",
    heights_equal = "row",
    bslib::card(
      style = 'background-color: #F7FAFD88;',
      bc_info,
      shiny::br(),
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
          ns("plot2"),
          title = "PC components",
          info.text = pcc.info,
          caption = pcc.info,
          options = NULL,
          height = c("100%","70vh")
        ),
        PlotModuleUI(
          ns("plot3"),
          title = "Covariate analysis",
          info.text = covariate.info,
          caption = covariate.info,
          options = covariate.options
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

upload_module_batchcorrect_server <- function(
      id,
      r_X,
      r_samples,
      r_contrasts,
      r_results,
      is.count = FALSE
  ) {
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
        X <- playbase::svdImpute2(X) ## standard impute
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
            selected = choices3[-1])
        } else {
          shiny::updateSelectInput(session, "tech_params", choices = "<none>",
            selected = "<none>")
        }
        
        if("statistical" %in% names(bc$params) && length(bc$params$statistical)) {
          bparams <- sort(bc$params$statistical)
          choices2 = c("<none>",bparams)
          shiny::updateSelectInput(session, "batch_params", choices = choices2,
            selected = choices2[-1])
        } else {
          shiny::updateSelectInput(session, "batch_params", choices = "<none>",
            selected = "<none>")
        }

        
        bc$choices <- list(
          batch_params = bparams,
          tech_params  = tparams
        )
        
        bc
      })
            
      output_plot1 <- function() {
        res <- r_results()
        shiny::req(res)
        X0 <- isolate(logX())
        samples <- isolate(r_samples())
        pheno <- res$pheno
        xlist <- res$xlist

        type = "umap"
        if(input$clust.plottype == "heatmap") {
          type = "heatmap"
          pos <- NULL
        } else {
          type = "umap"
          layout <- input$clust.plottype
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
          cex = 0.95
          )
      }
            
      output_plot2 <- function() {
        res <- r_results()          
        X0 <- isolate(logX())
        samples <- isolate(r_samples())
        xlist <- res$xlist

        kk <- colnames(xlist[[1]])
        X0 <- X0[,kk]
        samples <- samples[kk,]        
        
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

      output_plot3 <- function() {
        bc <- analyze_batch_effects()       
        shiny::req(bc)
        if(input$covariate.plottype == "Covariate plot") {
          par(mfrow=c(2,2), mar=c(3.3,4,2,2), mgp=c(2.2,1,0))
          playbase::bc.CovariateAnalysisPlot(bc, k=1:3, par=FALSE, col='#273871')
        } else {
          par(mfrow=c(1,1), mar=c(0,3,0,3))
          playbase::bc.plotCovariateHeatmap(bc)
        }
        
      }
        
      output_plot4 <- function() {
        res <- r_results()                    
        shiny::req(res)
        X0 <- isolate(logX())
        samples <- isolate(r_samples())
        pheno = NULL
        xlist <- res$xlist
        sel <- c("score","genes","gsets","SNR","pc1.ratio","silhouette")
        sel <- intersect(sel, colnames(res$scores)) 
        scores <- res$scores[,sel]
                
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
