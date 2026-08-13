## This file is part of the Omics Playground project.

qsee_outlier_server <- function(id, rX, rY) {
  shiny::moduleServer(id, function(input, output, session) {
    OmicsBoard("board", pgx = NULL, title = "Outlier analysis", infotext = NULL)
    is_visible <- board_is_visible(input, label = "qsee_outlier_server")
    observers <- board_observer_registry()
    observers$add(shiny::observeEvent(rY(), {
      shiny::updateSelectInput(session, "colorby", choices = colnames(rY()))
    }))
    get_result <- board_cache(
      is_visible, deps = function() list(rX(), rY()), label = "qsee_outlier_server",
      compute = function() {
        X <- rX(); Y <- rY(); shiny::req(X, Y)
        progress <- shiny::Progress$new(session, min = 0, max = 1); on.exit(progress$close())
        qsee_outlier_compute(X, Y, progress)
      }
    )

    render.outlier_zscores <- function() {
      res <- get_result(); ph <- input$colorby; shiny::req(ph); qsee_outlier_plot_zscores_plotly(res, ph)
    }
    render.outlier_pca <- function() {
      res <- get_result(); ph <- input$colorby; shiny::req(ph)
      qsee_outlier_plot_pca_plotly(res, ph, show_labels = isTRUE(input$show_labels))
    }

    PlotModuleServer(
      "outlier_zscores",
      plotlib = "plotly",
      func = render.outlier_zscores,
      add.watermark = FALSE
    )

    ## Heatmap uses custom renderer (iheatmapr objects from omicsplots::pgx.plot_heatmap
    ## do not play nicely with PlotModuleServer + iheatmapr::renderIheatmap which forces to_widget())
    qsee_plotly_hm_server(output, "outlier_heatmap", function() {
      res <- get_result(); shiny::req(res$heatX); qsee_outlier_plot_heatmap_plotly(res)
    })

    PlotModuleServer(
      "outlier_pca",
      plotlib = "plotly",
      func = render.outlier_pca,
      add.watermark = FALSE
    )

    board_pause_resume_observers(is_visible, observers, label = "qsee_outlier_server", start_paused = TRUE)
  })
}
