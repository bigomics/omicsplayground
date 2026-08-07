## This file is part of the Omics Playground project.

qsee_outlier_server <- function(id, rX, rY) {
  shiny::moduleServer(id, function(input, output, session) {
    is_visible <- qsee_is_visible(input, label = "qsee_outlier_server")
    shiny::observeEvent(rY(), {
      shiny::updateSelectInput(session, "colorby", choices = colnames(rY()))
    })
    get_result <- qsee_board_cache(
      is_visible, deps = function() list(rX(), rY()), label = "qsee_outlier_server",
      compute = function() {
        X <- rX(); Y <- rY(); shiny::req(X, Y)
        progress <- shiny::Progress$new(session, min = 0, max = 1); on.exit(progress$close())
        qsee_outlier_compute(X, Y, progress)
      }
    )
    output$outlier_zscores <- plotly::renderPlotly({
      res <- get_result(); ph <- input$colorby; shiny::req(ph); qsee_outlier_plot_zscores_plotly(res, ph)
    })
    qsee_plotly_hm_server(output, "outlier_heatmap", function() {
      res <- get_result(); shiny::req(res$heatX); qsee_outlier_plot_heatmap_plotly(res)
    })
    output$outlier_pca <- plotly::renderPlotly({
      res <- get_result(); ph <- input$colorby; shiny::req(ph); qsee_outlier_plot_pca_plotly(res, ph)
    })
  })
}
