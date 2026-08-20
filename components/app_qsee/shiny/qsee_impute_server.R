## This file is part of the Omics Playground project.

qsee_imputation_server <- function(id, rX, rY, purge = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    OmicsBoard("board", pgx = NULL, title = "Missing value analysis", infotext = NULL)
    is_visible <- bigdash::bd_is_visible(
      input, purge = qsee_resolve_purge(purge), label = "qsee_imputation_server"
    )
    redraw_tick <- bigdash::bd_redraw_tick(session = session)

    output$ui_output <- shiny::renderUI({
      qsee_imputation_ui_output(session$ns)
    })
    shiny::observeEvent(rY(), {
      shiny::updateSelectInput(session, "colorby", choices = colnames(rY()))
    })
    no_missing_warned <- shiny::reactiveVal(FALSE)
    shiny::observeEvent(is_visible(), {
      shiny::req(isTRUE(is_visible()))
      shiny::req(!no_missing_warned())
      rawX <- rX()
      shiny::req(rawX)
      if (!anyNA(rawX)) {
        shinyalert::shinyalert(
          text = "note, data has no missing values",
          type = "warning"
        )
        no_missing_warned(TRUE)
      }
    })
    tr_marlevel <- shiny::reactive(input$marlevel) %>% shiny::debounce(1000)
    tr_mnarlevel <- shiny::reactive(input$mnarlevel) %>% shiny::debounce(1000)
    ## Lazy: only the board's plot outputs read this, and Shiny suspends
    ## those while the tab is hidden. See qsee_visibility.R.
    get_result <- shiny::reactive({
      rawX <- rX(); shiny::req(rawX)
      message("[qsee_imputation_server] computing...")
      progress <- shiny::Progress$new(session, min = 0, max = 1); on.exit(progress$close())
      qsee_imputation_compute(rawX, tr_marlevel(), tr_mnarlevel(), progress)
    })

    render.pca_plots <- function() {
      res <- get_result(); Y <- rY(); ph <- input$colorby
      shiny::req(ph, ph %in% colnames(Y))
      qsee_imputation_plot_pca_plotly(res, Y, ph, show_labels = isTRUE(input$show_labels))
    }
    render.histograms <- function() {
      qsee_imputation_plot_histograms_plotly(get_result())
    }
    render.distributions <- function() {
      res <- get_result(); Y <- rY(); ph <- input$colorby
      shiny::req(ph, ph %in% colnames(Y)); qsee_imputation_plot_distributions_plotly(res, Y, ph)
    }
    render.validation_scatter <- function() {
      qsee_imputation_plot_validation_plotly(get_result())
    }

    heatmap_panels <- shiny::reactive({
      redraw_tick()
      qsee_imputation_plot_heatmaps_plotly(get_result())
    })
    qsee_plotly_hm_grid_server(output, "heatmap", heatmap_panels, n = 6L)

    PlotModuleServer(
      "pca_plots",
      plotlib = "plotly",
      func = render.pca_plots,
      add.watermark = FALSE
    )
    PlotModuleServer(
      "histograms",
      plotlib = "plotly",
      func = render.histograms,
      add.watermark = FALSE
    )
    PlotModuleServer(
      "distributions",
      plotlib = "plotly",
      func = render.distributions,
      add.watermark = FALSE
    )
    PlotModuleServer(
      "validation_scatter",
      plotlib = "plotly",
      func = render.validation_scatter,
      add.watermark = FALSE
    )
  })
}
