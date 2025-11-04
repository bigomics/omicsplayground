##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## =============================================================================
## ==================== BATCHCORRECT UI/SERVER =================================
## =============================================================================

applet_batchcorrect_inputs <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tagList(
    shiny::radioButtons(ns("clust.plottype"), "Type;", c("umap", "heatmap")),
    shiny::radioButtons( ns("covariate.plottype"), "Type;",
      c("Covariate plot", "Correlation heatmap")),
    shiny::selectizeInput(ns("tech_params"), "Technical covariates:",
      choices = "<none>", multiple = TRUE ),
    shiny::selectizeInput(ns("batch_params"), "Batch covariates:",
      choices = "<none>", multiple = TRUE),
    shiny::actionButton(ns("recompute_button"), "Recompute",
      class = "btn-sm btn-primary mt-3", width = "100%"
    )
  )
  
}

applet_batchcorrect_ui <- function(id, as="ui") {
  ns <- shiny::NS(id)

  clust.infotext <-
    "Clustering of samples before ('uncorrected') and after the different batch-effect correction methods. After batch-effect correction clustering should improve. The silhouette score gives an indication of the clustering performance of the method."

  pcc.info <- "PC analysis by covariate (class). The heights of the bars correspond to the relative contribution of that covariate to the three PC, as measured by an F-test. "

  covariate.info <-
    "Analysis of variables by plotting their significance in correlation with the phenotype against their significance in correlation with a principal component (PC) vector. Strong model variables are situate 'top right'. Batch effect variables with high PC correlation but low phenotype correlation are on the 'top left'. A well-designed experiment shows strong model variables in PC1, else it may be a sign of significant batch-effects."

  modules <- tagList(
    PlotModuleUI(
      ns("plot1"),
      title = "Clustering",
      info.text = clust.infotext,
      caption = clust.infotext,
      #options = clust.options,
      height = c("100%", "70vh")
    ),
    PlotModuleUI(
      ns("plot2"),
      title = "PC components",
      info.text = pcc.info,
      caption = pcc.info,
      options = NULL,
      height = c("100%", "70vh")
    ),
    PlotModuleUI(
      ns("plot3"),
      title = "Covariate analysis",
      info.text = covariate.info,
      caption = covariate.info,
      #options = covariate.options
    ),
    PlotModuleUI(
      ns("plot4"),
      title = "Statistics and score",
      options = NULL
    )
  )

  modules_info <- list(
    "plot1" = clust.infotext,
    "plot2" = pcc.info,
    "plot3" = covariate.info,
    "plot4" = ""
  )

  if(as=="taglist") return(modules)
  if(as=="info") return(modules_info)

  bslib::layout_columns(
    col_widths = 6,
    height = "calc(100vh - 200px)",
    heights_equal = "row",
    gap = "35px",
    !!!modules
  )
  
}


applet_batchcorrect_server <- function(id, pgx) {
  shiny::moduleServer(
    id, function(input, output, session) {

      analyze_batch_effects <- reactive({

        bc <- playbase::detectBatchEffects(
          X = pgx$X,
          samples = pgx$samples,
          pheno = NULL,
          contrasts = pgx$contrasts,
          ## params = c("statistical","technical","pca"),
          params = c("statistical", "technical"),
          p.pca = 0.5,
          p.pheno = 0.05,
          k.pca = 10,
          nv = 2,
          xrank = NULL
        )

        bparams <- c()
        tparams <- c()

        if ("technical" %in% names(bc$params) && length(bc$params$technical)) {
          tparams <- sort(unique(sub("[.].*", "", bc$params$technical)))
          choices3 <- c("<none>", paste0("<", tparams, ">"))
          shiny::updateSelectInput(session, "tech_params",
            choices = choices3,
            selected = choices3[-1]
          )
        } else {
          shiny::updateSelectInput(session, "tech_params",
            choices = "<none>",
            selected = "<none>"
          )
        }

        if ("statistical" %in% names(bc$params) && length(bc$params$statistical)) {
          bparams <- sort(bc$params$statistical)
          choices2 <- c("<none>", bparams)
          shiny::updateSelectInput(session, "batch_params",
            choices = choices2,
            selected = choices2[-1]
          )
        } else {
          shiny::updateSelectInput(session, "batch_params",
            choices = "<none>",
            selected = "<none>"
          )
        }

        bc$choices <- list(
          batch_params = bparams,
          tech_params  = tparams
        )

        bc

      })
      
      compare_bc_methods <- reactive({

        pheno = pgx$samples[,"time"]
        batch.pars = NULL
        pars <- playbase::get_model_parameters(
          pgx$X, pgx$samples, pheno = pheno, contrasts = NULL)
    
        res <- playbase::compare_batchcorrection_methods(
          pgx$X,
          pgx$samples,
          pheno = pars$pheno,
          contrasts = NULL,
          batch.pars = pars$batch.pars,
          clust.method = "pca",
          evaluate = TRUE, ## no score computation
          #xlist.init = xlist.init,
          ntop = 1000
        )
        
        res
      })

      
      output_plot1 <- function() {

        res <- compare_bc_methods()
        pheno <- res$pheno
        xlist <- res$xlist

        pos <- res$pos
        type <- "umap"
        type <- input$clust.plottype
        
        playbase::bc.plotResults(
          X = pgx$X,
          xlist = xlist,
          pos = pos,
          pheno = pheno,
          samples = pgx$samples,
          type = type,
          scores = NULL,
          nmax = 1000,
          cex = 0.95
        )

      }

      output_plot2 <- function() {

        res <- compare_bc_methods()
        pheno <- res$pheno
        xlist <- res$xlist

        kk <- colnames(xlist[[1]])
        X1 <- pgx$X[, kk]
        samples1 <- pgx$samples[kk, ]

        playbase::bc.plotResults(
          X = X1,
          xlist = xlist,
          pos = NULL,
          pheno = NULL,
          samples = samples1,
          type = "pc",
          scores = NULL,
          nmax = 1000
        )

      }

      output_plot3 <- function() {
        bc <- analyze_batch_effects()
        shiny::req(bc)
        if(input$covariate.plottype == "Covariate plot") {      
          par(mfrow = c(2, 2), mar = c(3.3, 4, 2, 2), mgp = c(2.2, 1, 0))
          playbase::bc.CovariateAnalysisPlot(bc, k = 1:3, par=FALSE, col="#273871")
        } else {
          par(mfrow = c(1, 1), mar = c(0, 3, 0, 3))
          playbase::bc.plotCovariateHeatmap(bc)
        }
      }

      output_plot4 <- function() {

        res <- compare_bc_methods()
        shiny::req(res)

        xlist <- res$xlist
        sel <- c("score", "genes", "gsets", "SNR", "pc1.ratio", "silhouette")
        sel <- intersect(sel, colnames(res$scores))
        scores <- res$scores[, sel]

        playbase::bc.plotResults(
          X = pgx$X,
          xlist = xlist,
          pos = NULL,
          pheno = NULL,
          samples = pgx$samples,
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
        res = c(65, 110),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      ##return(correctedX) ## pointing to reactive
    } ## end-of-server
  )
}



applet_batchcorrect <- list(
  title = "Batchcorrect Applet",
  inputs = applet_batchcorrect_inputs,
  ui = applet_batchcorrect_ui,  
  server = applet_batchcorrect_server
)
