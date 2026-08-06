## This file is part of the Omics Playground project.

qsee_imputation_server <- function(id, rX, rY) {
  shiny::moduleServer(id, function(input, output, session) {
    is_visible <- qsee_is_visible(input, label = "qsee_imputation_server")
    shiny::observeEvent(rY(), {
      shiny::updateSelectInput(session, "colorby", choices = colnames(rY()))
    })
    tr_marlevel <- shiny::reactive(input$marlevel) %>% shiny::debounce(1000)
    tr_mnarlevel <- shiny::reactive(input$mnarlevel) %>% shiny::debounce(1000)
    get_result <- qsee_board_cache(
      is_visible,
      deps = function() list(rX(), tr_marlevel(), tr_mnarlevel()),
      label = "qsee_imputation_server",
      compute = function() {
        rawX <- rX(); shiny::req(rawX)
        progress <- shiny::Progress$new(session, min = 0, max = 1); on.exit(progress$close())
        qsee_imputation_compute(rawX, tr_marlevel(), tr_mnarlevel(), progress)
      }
    )
    output$pca_plots <- shiny::renderPlot({
      res <- get_result(); Y <- rY(); ph <- input$colorby
      shiny::req(ph, ph %in% colnames(Y)); qsee_imputation_plot_pca(res, Y, ph)
    })
    output$heatmap_plots <- shiny::renderPlot({ qsee_imputation_plot_heatmaps(get_result()) })
    output$histograms <- shiny::renderPlot({ qsee_imputation_plot_histograms(get_result()) })
    output$distributions <- shiny::renderPlot({
      res <- get_result(); Y <- rY(); ph <- input$colorby
      shiny::req(ph, ph %in% colnames(Y)); qsee_imputation_plot_distributions(res, Y, ph)
    })
    output$validation_scatter <- shiny::renderPlot({ qsee_imputation_plot_validation(get_result()) })
  })
}
