## This file is part of the Omics Playground project.

qsee_normalization_server <- function(id, rX, rY) {
  shiny::moduleServer(id, function(input, output, session) {

    OmicsBoard("board", pgx=NULL, title="Normalization", infotext = NULL) 

    is_visible <- qsee_is_visible(input, label = "qsee_normalization_server")
    shiny::observeEvent(rY(), {
      shiny::updateSelectInput(session, "colorby", choices = colnames(rY()))
    })

    get_rawX <- shiny::reactive({
      rawX <- rX(); shiny::req(rawX)
      amount <- input$noise
      if (is.null(amount)) amount <- 1
      qsee_normalization_add_noise(rawX, amount)
    })

    get_result <- qsee_board_cache(
      is_visible, deps = function() get_rawX(), label = "qsee_normalization_server",
      compute = function() {
        rawX <- get_rawX(); shiny::req(rawX)
        progress <- shiny::Progress$new(session, min = 0, max = 1); on.exit(progress$close())
        progress$set(message = "Compute normalization...", value = 0.33)
        qsee_normalization_compute(rawX, progress)
      }
    )

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
