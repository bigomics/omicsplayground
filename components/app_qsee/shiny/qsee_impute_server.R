## This file is part of the Omics Playground project.

qsee_imputation_server <- function(id, rX, rY) {
  shiny::moduleServer(id, function(input, output, session) {
    OmicsBoard("board", pgx = NULL, title = "Missing value analysis", infotext = NULL)
    is_visible <- board_is_visible(input, label = "qsee_imputation_server")
    observers <- board_observer_registry()
    observers$add(shiny::observeEvent(rY(), {
      shiny::updateSelectInput(session, "colorby", choices = colnames(rY()))
    }))
    tr_marlevel <- shiny::reactive(input$marlevel) %>% shiny::debounce(1000)
    tr_mnarlevel <- shiny::reactive(input$mnarlevel) %>% shiny::debounce(1000)
    get_result <- board_cache(
      is_visible,
      deps = function() list(rX(), tr_marlevel(), tr_mnarlevel()),
      label = "qsee_imputation_server",
      compute = function() {
        rawX <- rX(); shiny::req(rawX)
        progress <- shiny::Progress$new(session, min = 0, max = 1); on.exit(progress$close())
        qsee_imputation_compute(rawX, tr_marlevel(), tr_mnarlevel(), progress)
      }
    )

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

    board_pause_resume_observers(is_visible, observers, label = "qsee_imputation_server")
  })
}
