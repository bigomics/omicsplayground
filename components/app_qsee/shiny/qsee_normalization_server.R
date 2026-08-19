## This file is part of the Omics Playground project.

qsee_normalization_server <- function(id, rX, rY, purge = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    OmicsBoard("board", pgx=NULL, title="Normalization", infotext = NULL) 

    is_visible <- bigdash::bd_is_visible(
      input, purge = qsee_resolve_purge(purge), label = "qsee_normalization_server"
    )
    ## Shiny suspends this output while the board's tab is hidden, so the
    ## body is only built on the first visit and then kept in the DOM.
    output$ui_output <- shiny::renderUI({
      qsee_normalization_ui_output(session$ns)
    })
    shiny::observeEvent(rY(), {
      shiny::updateSelectInput(session, "colorby", choices = colnames(rY()))
    })

    get_rawX <- shiny::reactive({
      rawX <- rX(); shiny::req(rawX)
      amount <- input$noise
      if (is.null(amount)) amount <- 1
      qsee_normalization_add_noise(rawX, amount)
    })

    ## Lazy: only the board's plot outputs read this, and Shiny suspends
    ## those while the tab is hidden. See qsee_visibility.R.
    get_result <- shiny::reactive({
      rawX <- get_rawX(); shiny::req(rawX)
      message("[qsee_normalization_server] computing...")
      progress <- shiny::Progress$new(session, min = 0, max = 1); on.exit(progress$close())
      progress$set(message = "Compute normalization...", value = 0.33)
      qsee_normalization_compute(rawX, progress)
    })

    render.box_plots <- function() {
      res <- get_result(); Y <- rY(); ph <- input$colorby
      shiny::req(ph, ph %in% colnames(Y)); qsee_normalization_plot_boxplots_plotly(res, Y, ph)
    }
    render.histograms <- function() {
      res <- get_result(); Y <- rY(); ph <- input$colorby
      shiny::req(ph, ph %in% colnames(Y)); qsee_normalization_plot_histograms_plotly(res, Y, ph)
    }
    render.pca_plots <- function() {
      res <- get_result(); Y <- rY(); ph <- input$colorby
      shiny::req(ph, ph %in% colnames(Y))
      qsee_normalization_plot_pca_plotly(res, Y, ph, show_labels = isTRUE(input$show_labels))
    }

    PlotModuleServer(
      "box_plots",
      plotlib = "plotly",
      func = render.box_plots,
      add.watermark = FALSE
    )
    PlotModuleServer(
      "histograms",
      plotlib = "plotly",
      func = render.histograms,
      add.watermark = FALSE
    )
    PlotModuleServer(
      "pca_plots",
      plotlib = "plotly",
      func = render.pca_plots,
      add.watermark = FALSE
    )
  })
}
