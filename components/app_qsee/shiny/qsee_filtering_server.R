## This file is part of the Omics Playground project.

qsee_filtering_server <- function(id, rX, rY) {
  shiny::moduleServer(id, function(input, output, session) {
    OmicsBoard("board", pgx = NULL, title = "SD filtering", infotext = NULL)
    is_visible <- bigdash::bd_is_visible(
      input, purge = qsee_purge_enabled(), label = "qsee_filtering_server"
    )

    output$ui_output <- shiny::renderUI({
      qsee_filtering_ui_output(session$ns)
    })

    shiny::observeEvent(rY(), {
      cols <- colnames(rY())
      shiny::updateSelectInput(session, "colorby", choices = cols, selected = cols[1])
    }, ignoreNULL = TRUE)

    ## A plain observe(), not observeEvent(rX(), ...): there the event
    ## expression *is* the pull, so a req() in the handler comes too late.
    ## Here req() aborts before rX() is read, and the observer re-runs on the
    ## first visit because it still depends on is_visible().
    shiny::observe({
      shiny::req(is_visible())
      X <- rX()
      shiny::req(X)
      nmax <- max(1L, nrow(X))
      shiny::updateSliderInput(session, "threshold", min = 1, max = nmax, value = min(1000L, nmax), step = 1)
    })

    get_results <- shiny::reactive({
      X <- rX(); Y <- rY(); shiny::req(X, Y)
      message("[qsee_filtering_server] computing...")
      qsee_filtering_compute(X, Y)
    })

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
  })
}
