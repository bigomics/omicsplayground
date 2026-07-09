##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Was: board.upload/upload_module_batchcorrect.R


qsee_server <- function(id,
                        r_X,
                        r_samples,
                        r_contrasts,
                        r_results,   ## BATCH CORRECT results???
                        is.count = FALSE) {
  shiny::moduleServer(
    id,
    function(input, output, session) {

      ## also return object
      correctedX <- reactiveVal(NULL)
      observeEvent(r_X(), correctedX(r_X()))

      logX <- eventReactive(
        {
          list(r_X())
        },
        {
          shiny::validate(shiny::need(!is.null(r_X()), "no data. please upload."))
          shiny::validate(shiny::need(!is.null(nrow(r_X())), "no data. please upload."))
          X <- playbase::counts.mergeDuplicateFeatures(r_X())
          X <- limma::normalizeQuantiles(playbase::logCPM(X)) ## standard normalization
          X <- playbase::svdImpute2(X) ## standard impute
          X
        }
      )

      analyze_batch_effects <- reactive({
        ## shiny::req(r_contrasts(), r_samples(), r_X())
        shiny::req(dim(r_contrasts()), dim(r_samples()), dim(logX()))
        samples <- r_samples()
        contrasts <- r_contrasts()
        X <- logX()

        shiny::validate(shiny::need(
          !is.null(r_X()) && nrow(r_X()),
          "No counts data. Please upload."
        ))
        shiny::validate(shiny::need(
          !is.null(r_contrasts()) && nrow(r_contrasts()),
          "No contrasts defined. Please create contrasts."
        ))

        bc <- playbase::detectBatchEffects(
          X = X,
          samples = samples,
          pheno = NULL,
          contrasts = contrasts,
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

      output_plot1 <- function() {
        res <- r_results()
        shiny::req(res)
        X0 <- isolate(logX())
        samples <- isolate(r_samples())
        pheno <- res$pheno
        xlist <- res$xlist

        type <- "umap"
        if (input$clust.plottype == "heatmap") {
          type <- "heatmap"
          pos <- NULL
        } else {
          type <- "umap"
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
        X0 <- X0[, kk]
        samples <- samples[kk, ]

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
        if (input$covariate.plottype == "Covariate plot") {
          par(mfrow = c(2, 2), mar = c(3.3, 4, 2, 2), mgp = c(2.2, 1, 0))
          playbase::bc.CovariateAnalysisPlot(bc, k = 1:3, par = FALSE, col = "#273871")
        } else {
          par(mfrow = c(1, 1), mar = c(0, 3, 0, 3))
          playbase::bc.plotCovariateHeatmap(bc)
        }
      }

      output_plot4 <- function() {
        res <- r_results()
        shiny::req(res)
        X0 <- isolate(logX())
        samples <- isolate(r_samples())
        pheno <- NULL
        xlist <- res$xlist
        sel <- c("score", "genes", "gsets", "SNR", "pc1.ratio", "silhouette")
        sel <- intersect(sel, colnames(res$scores))
        scores <- res$scores[, sel]

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


      return(correctedX) ## pointing to reactive
    } ## end-of-server
  )
}
