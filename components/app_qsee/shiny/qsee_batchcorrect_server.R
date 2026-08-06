##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

qsee_bsee_server <- function(id, rX, rY) {
  shiny::moduleServer(
    id,
    function(input, output, session) {

      ## `input$is_visible` is reported by qsee_visibility_probe() in the UI.
      is_visible <- qsee_is_visible(input, label = "qsee_bsee_server")

      observeEvent( list(rY()), {
        Y <- rY()
        shiny::req(Y)
        ny <- apply(Y, 2, function(y) max(table(y)))
        sel.main <- which.max( ny * ( ny < nrow(Y) & ny > 1 ) )
        shiny::updateSelectInput(session, "main_param",
          choices = colnames(Y),
          selected = colnames(Y)[sel.main]
        )
      })
      
      ## Main precomputation. Visibility gates *when* work may run; cache
      ## invalidates only when inputs / Recompute change — not on tab
      ## leave/return.
      get_results <- qsee_board_cache(
        is_visible,
        deps = function() list(rX(), rY(), input$recompute_button),
        label = "qsee_bsee_server",
        compute = function() {
          X <- rX()
          samples <- rY()
          shiny::req(X, samples)

          progress <- shiny::Progress$new(session, min = 0, max = 1)
          on.exit(progress$close())
          progress$set(message = paste("Normalizing..."), value = 0.1)

          pheno <- colnames(samples)[1]  ## NEED UPDATE!!!
          pheno <- input$main_param
          shiny::req(pheno)
          
          res <- bsee_compute_batchcorrect(X, samples, pheno, progress = progress) 
          return(res)
        }
      )

      ## Update selectInputs when results are available and the tab is visible.
      ## Because get_results() now gates on is_visible() internally, we no
      ## longer need the req() here for performance reasons (but we keep
      ## it for clarity).
      observe({
        req(is_visible())
        res <- get_results()

        bparams <- c()
        tparams <- c()
        bc <- res$effects
                
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
        
        methods <- c("ComBat", "limma", "RUV", "SVA", "NPM")
        shiny::updateSelectInput(
          session,
          "bec_method",
          choices = methods
        )

        pp <- colnames(res$samples)
        shiny::updateSelectInput( session, "clust.colorby", choices = pp)
        
      })
           
      render.plot_pca_vs_methods <- function() {
        res <- get_results()
        ##pheno.var <- input$clust.colorby
        pheno.var <- input$main_param
        shiny::req(res, pheno.var)
        bsee.plot_pca_vs_methods(res, pheno.var, cex=1.2) 
      }

      render.plot_heatmap_vs_methods <- function() {
        res <- get_results()
        shiny::req(res)
        bsee.plot_heatmap_vs_methods(res) 
      }

      render.plot_covariate_correlation_heatmap <- function() {
        res <- get_results()
        shiny::req(res)
        bsee.plot_covariate_correlation_heatmap(res) 
      }

      render.plot_covariate_analysis <- function() {
        res <- get_results()
        shiny::req(res)
        bsee.plot_covariate_analysis(res) 
      }

      render.plot_scores <- function() {
        res <- get_results()
        shiny::req(res)
        bsee.plot_scores(res) 
      }

      render.plot_pvca_by_phenotype <- function() {
        res <- get_results()
        shiny::req(res)
        bsee.plot_pvca_by_phenotype(res) 
      }

      render.plot_pvca_by_component <- function() {
        res <- get_results()
        shiny::req(res)
        bsee.plot_pvca_by_component(res) 
      }
      
      PlotModuleServer(
        "plot1",
        plotlib = "base",
        func = render.plot_pca_vs_methods,
        res = c(70, 110),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot2",
        plotlib = "base",
        func = render.plot_heatmap_vs_methods,
        res = c(65, 110),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot3",
        plotlib = "base",
        func = render.plot_scores,
        res = c(65, 110),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot4",
        plotlib = "base",
        func = render.plot_covariate_correlation_heatmap,
        res = c(70, 110),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      ##-------------------------------
      
      PlotModuleServer(
        "plot5",
        plotlib = "base",
        func = render.plot_pvca_by_phenotype,
        res = c(65, 110),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot6",
        plotlib = "base",
        func = render.plot_pvca_by_component,
        res = c(65, 110),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot7",
        plotlib = "base",
        func = render.plot_covariate_analysis,
        res = c(65, 110),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )
          
      return(NULL) ## pointing to reactive
    } ## end-of-server
  )
}
