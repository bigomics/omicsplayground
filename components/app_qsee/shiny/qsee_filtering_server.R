## This file is part of the Omics Playground project.

qsee_filtering_server <- function(id, rX, rY) {
  shiny::moduleServer(id, function(input, output, session) {
    OmicsBoard("board", pgx = NULL, title = "SD filtering", infotext = NULL)
    is_visible <- board_is_visible(input, label = "qsee_filtering_server")
    observers <- board_observer_registry()
    observers$add(shiny::observeEvent(rY(), {
      cols <- colnames(rY())
      shiny::updateSelectInput(session, "colorby", choices = cols, selected = cols[1])
    }, ignoreNULL = TRUE))
    observers$add(shiny::observeEvent(rX(), {
      nmax <- max(1L, nrow(rX()))
      shiny::updateSliderInput(session, "threshold", min = 1, max = nmax, value = min(1000L, nmax), step = 1)
    }, ignoreNULL = TRUE))
    get_results <- board_cache(
      is_visible, deps = function() list(rX(), rY()), label = "qsee_filtering_server",
      compute = function() { X <- rX(); Y <- rY(); shiny::req(X, Y); qsee_filtering_compute(X, Y) }
    )

    render.pca_vs_topsd <- function() {
      res <- get_results(); ph <- input$colorby
      shiny::req(ph)
      qsee_filtering_plot_pca_vs_topsd_plotly(res, ph, show_labels = isTRUE(input$show_labels))
    }
    render.variance_vs_topsd <- function() {
      res <- get_results(); top <- input$threshold
      shiny::req(top); qsee_filtering_plot_variance_vs_topsd_plotly(res, top)
    }
    render.sd_histogram <- function() {
      res <- get_results(); top <- input$threshold
      shiny::req(top); qsee_filtering_plot_sd_histogram_plotly(res, top)
    }

    PlotModuleServer(
      "pca_vs_topsd",
      plotlib = "plotly",
      func = render.pca_vs_topsd,
      add.watermark = FALSE
    )
    PlotModuleServer(
      "variance_vs_topsd",
      plotlib = "plotly",
      func = render.variance_vs_topsd,
      add.watermark = FALSE
    )
    PlotModuleServer(
      "sd_histogram",
      plotlib = "plotly",
      func = render.sd_histogram,
      add.watermark = FALSE
    )

    board_pause_resume_observers(is_visible, observers, label = "qsee_filtering_server")
  })
}
