## This file is part of the Omics Playground project.

qsee_normalization_server <- function(id, rX, rY) {
  shiny::moduleServer(id, function(input, output, session) {
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

    output$box_plots <- plotly::renderPlotly({
      res <- get_result(); Y <- rY(); ph <- input$colorby
      shiny::req(ph, ph %in% colnames(Y)); qsee_normalization_plot_boxplots_plotly(res, Y, ph)
    })
    output$histograms <- plotly::renderPlotly({
      res <- get_result(); Y <- rY(); ph <- input$colorby
      shiny::req(ph, ph %in% colnames(Y)); qsee_normalization_plot_histograms_plotly(res, Y, ph)
    })
    output$pca_plots <- plotly::renderPlotly({
      res <- get_result(); Y <- rY(); ph <- input$colorby
      shiny::req(ph, ph %in% colnames(Y)); qsee_normalization_plot_pca_plotly(res, Y, ph)
    })
  })
}
