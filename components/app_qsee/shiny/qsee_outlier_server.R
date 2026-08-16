## This file is part of the Omics Playground project.

qsee_outlier_server <- function(id, rX, rY) {
  shiny::moduleServer(id, function(input, output, session) {
    OmicsBoard("board", pgx = NULL, title = "Outlier analysis", infotext = NULL)
    is_visible <- qsee_is_visible(input, label = "qsee_outlier_server")
    redraw_tick <- qsee_plotly_purge(is_visible, session, label = "qsee_outlier_server")

    output$ui_output <- shiny::renderUI({
      qsee_outlier_ui_output(session$ns)
    })

    shiny::observeEvent(rY(), {
      shiny::updateSelectInput(session, "colorby", choices = colnames(rY()))
    })

    ## Lazy: only the board's plot outputs read this, and Shiny suspends
    ## those while the tab is hidden. See qsee_visibility.R.
    get_result <- shiny::reactive({
      X <- rX(); Y <- rY(); shiny::req(X, Y)
      message("[qsee_outlier_server] computing...")
      progress <- shiny::Progress$new(session, min = 0, max = 1); on.exit(progress$close())
      qsee_outlier_compute(X, Y, progress)
    })

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
      func = qsee_with_redraw(redraw_tick, render.outlier_zscores),
      add.watermark = FALSE
    )

    ## Heatmap uses custom renderer (iheatmapr objects from omicsplots::pgx.plot_heatmap
    ## do not play nicely with PlotModuleServer + iheatmapr::renderIheatmap which forces to_widget())
    qsee_plotly_hm_server(output, "outlier_heatmap", function() {
      redraw_tick()
      res <- get_result(); shiny::req(res$heatX); qsee_outlier_plot_heatmap_plotly(res)
    })

    PlotModuleServer(
      "outlier_pca",
      plotlib = "plotly",
      func = qsee_with_redraw(redraw_tick, render.outlier_pca),
      add.watermark = FALSE
    )
  })
}
