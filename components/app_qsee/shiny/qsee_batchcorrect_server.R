##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

qsee_bsee_server <- function(id, rX, rY) {
  shiny::moduleServer(
    id,
    function(input, output, session) {

      OmicsBoard("board", pgx = NULL, title = "Batch-effects", infotext = NULL)
      ## Visibility is reported by bigdash::bd_visibility_probe() in the UI.
      is_visible <- bigdash::bd_is_visible(
        input, purge = qsee_purge_enabled(), label = "qsee_bsee_server"
      )
      redraw_tick <- bigdash::bd_redraw_tick(session = session)
      ## Shiny suspends this output while the board's tab is hidden, so the
      ## body is only built on the first visit and then kept in the DOM.
      output$ui_output <- shiny::renderUI({
        qsee_bsee_ui_output(session$ns)
      })

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
      
      ## Main precomputation. eventReactive, not reactive: the body reads
      ## input$main_param without depending on it, so changing the main
      ## parameter does not kick off a batch correction -- only the inputs
      ## listed below do, i.e. new data or the Recompute button.
      ##
      ## Lazy, like any reactive: the plot outputs that read it are suspended
      ## while the board is hidden, so nothing runs until the tab is opened
      ## (but see the observe() below).
      get_results <- shiny::eventReactive(
        list(rX(), rY(), input$recompute_button),
        {
          X <- rX()
          samples <- rY()
          shiny::req(X, samples)

          message("[qsee_bsee_server] computing...")
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
      ##
      ## DO NOT REMOVE the req(is_visible()) below. This is the only eager
      ## consumer of get_results() -- every other reader is a plot output,
      ## which Shiny suspends while the board is hidden. Without the guard
      ## this observer pulls the cache as soon as rX()/rY() arrive and the
      ## whole batch correction runs even if the tab is never opened.
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
        redraw_tick()
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
    } ## end-of-server
  )
}
