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
    output$outlier_zscores <- shiny::renderPlot({
      res <- get_result(); ph <- input$colorby; shiny::req(ph); qsee_outlier_plot_zscores(res, ph)
    })
    output$outlier_heatmap <- shiny::renderPlot({
      res <- get_result(); shiny::req(res$heatX); qsee_outlier_plot_heatmap(res)
    })
    output$outlier_pca <- shiny::renderPlot({
      res <- get_result(); ph <- input$colorby; shiny::req(ph); qsee_outlier_plot_pca(res, ph)
    })
  })
}
