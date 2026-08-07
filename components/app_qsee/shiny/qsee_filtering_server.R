## This file is part of the Omics Playground project.

qsee_filtering_server <- function(id, rX, rY) {
  shiny::moduleServer(id, function(input, output, session) {
    is_visible <- qsee_is_visible(input, label = "qsee_filtering_server")
    shiny::observeEvent(rY(), {
      cols <- colnames(rY())
      shiny::updateSelectInput(session, "colorby", choices = cols, selected = cols[1])
    }, ignoreNULL = TRUE)
    shiny::observeEvent(rX(), {
      nmax <- max(1L, nrow(rX()))
      shiny::updateSliderInput(session, "threshold", min = 1, max = nmax, value = min(1000L, nmax), step = 1)
    }, ignoreNULL = TRUE)
    get_results <- qsee_board_cache(
      is_visible, deps = function() list(rX(), rY()), label = "qsee_filtering_server",
      compute = function() { X <- rX(); Y <- rY(); shiny::req(X, Y); qsee_filtering_compute(X, Y) }
    )
    output$pca_vs_topsd <- plotly::renderPlotly({
      res <- get_results(); ph <- input$colorby; shiny::req(ph); qsee_filtering_plot_pca_vs_topsd_plotly(res, ph)
    })
    output$variance_vs_topsd <- plotly::renderPlotly({
      res <- get_results(); top <- input$threshold; shiny::req(top); qsee_filtering_plot_variance_vs_topsd_plotly(res, top)
    })
    output$sd_histogram <- plotly::renderPlotly({
      res <- get_results(); top <- input$threshold; shiny::req(top); qsee_filtering_plot_sd_histogram_plotly(res, top)
    })
  })
}
