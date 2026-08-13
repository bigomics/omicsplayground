##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

qsee_bsee_server <- function(id, rX, rY) {
  shiny::moduleServer(
    id,
    function(input, output, session) {

      OmicsBoard("board", pgx = NULL, title = "Batch-effects", infotext = NULL)
      ## `input$is_visible` is reported by board_visibility_probe() in the UI.
      is_visible <- board_is_visible(input, label = "qsee_bsee_server")
      observers <- board_observer_registry()

      ## NOT suspended: compute() below reads input$main_param via an
      ## isolate()'d call inside board_cache, outside of deps() tracking.
      ## If this setter were suspended while hidden, the first visibility
      ## flip would resume it too late -- the (never-suspended) compute
      ## observer reacts to is_visible() directly and would already have run
      ## with an empty main_param, req()-failing silently until a *second*
      ## visibility toggle gave it another chance. Must stay always-live.
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
      get_results <- board_cache(
        is_visible,
        deps = function() list(rX(), rY(), input$recompute_button, input$main_param),
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

      ## Update selectInputs when results are available. Suspended/resumed
      ## explicitly by board_pause_resume_observers() below, so no req() on
      ## is_visible() is needed here -- it simply won't run while hidden.
      observers$add(observe({
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

      }))

      render.plot_pca_vs_methods <- function() {
        res <- get_results()
        pheno.var <- input$main_param
        shiny::req(res, pheno.var)
        bsee.plot_pca_vs_methods_plotly(res, pheno.var, show_labels = isTRUE(input$show_labels))
      }

      render.plot_covariate_correlation_heatmap <- function() {
        res <- get_results()
        shiny::req(res)
        bsee.plot_covariate_correlation_heatmap_plotly(res)
      }

      render.plot_covariate_analysis <- function() {
        res <- get_results()
        shiny::req(res)
        bsee.plot_covariate_analysis_plotly(res)
      }

      render.plot_scores <- function() {
        res <- get_results()
        shiny::req(res)
        bsee.plot_scores_plotly(res)
      }

      render.plot_pvca_by_phenotype <- function() {
        res <- get_results()
        shiny::req(res)
        bsee.plot_pvca_by_phenotype_plotly(res)
      }

      render.plot_pvca_by_component <- function() {
        res <- get_results()
        shiny::req(res)
        bsee.plot_pvca_by_component_plotly(res)
      }

      PlotModuleServer(
        "plot1",
        plotlib = "plotly",
        func = render.plot_pca_vs_methods,
        add.watermark = FALSE
      )

      heatmap_panels <- shiny::reactive({
        res <- get_results()
        shiny::req(res)
        bsee.plot_heatmap_vs_methods_plotly(res)
      })
      qsee_plotly_hm_grid_server(output, "bsee", heatmap_panels, n = 6L)

      PlotModuleServer(
        "plot3",
        plotlib = "plotly",
        func = render.plot_scores,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot4",
        plotlib = "plotly",
        func = render.plot_covariate_correlation_heatmap,
        add.watermark = FALSE
      )

      ##-------------------------------

      PlotModuleServer(
        "plot5",
        plotlib = "plotly",
        func = render.plot_pvca_by_phenotype,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot6",
        plotlib = "plotly",
        func = render.plot_pvca_by_component,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot7",
        plotlib = "plotly",
        func = render.plot_covariate_analysis,
        add.watermark = FALSE
      )

      board_pause_resume_observers(is_visible, observers, label = "qsee_bsee_server", start_paused = TRUE)
    } ## end-of-server
  )
}
