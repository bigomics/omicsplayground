PCAexplorer_biplot_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_columns(
    col_widths = c(7, 5),
    bslib::navset_card_tab(
      full_screen = TRUE,
      title = "Biplot",
      bslib::nav_spacer(),
      bslib::nav_panel("base", plotOutput(ns("biplot"), height = "800px")),
      bslib::nav_panel("plotly", plotly::plotlyOutput(ns("biplot_2d"), height = "800px")),
      bslib::nav_panel("3D", plotly::plotlyOutput(ns("biplot_3d"), height = "800px")),
      bslib::nav_panel("table", dataTableOutput(ns("biplot_table")))
    ),
    bslib::layout_columns(
      col_widths = 12,
      row_heights = c(1, 1),
      bslib::card(full_screen = TRUE, plotOutput(ns("variance_proportion"), height = "400px")),
      bslib::card(full_screen = TRUE, plotOutput(ns("variance_cumulative"), height = "400px"))
    )
  )
}

PCAexplorer_biplot_module <- function(id, get_pcaX,
                                      colorby, numarrows, show_arrows,
                                      arrowtype, ellipse,
                                      var_scale = shiny::reactive(0.65)) {
  shiny::moduleServer(
    id,
    function(input, output, session) {

      ## 0 arrows when show_arrows is off
      n_arrows <- shiny::reactive({
        if (isTRUE(show_arrows())) as.integer(numarrows()) else 0L
      })

      vscale <- shiny::reactive({
        v <- suppressWarnings(as.numeric(var_scale()))
        if (!is.finite(v) || v <= 0) 0.65 else v
      })

      output$biplot_table <- renderDataTable({
        res <- get_pcaX()
        dt <- round(res$v, 4)

        DT::datatable(
          data = dt,
          class = "compact hover",
          rownames = TRUE,
          extensions = c("Buttons", "Scroller"),
          plugins = "scrollResize",
          selection = list(
            mode = "single",
            target = "row",
            selected = 1
          ),
          options = list(
            dom = "lfrtipB",
            buttons = c('copy', 'csv', 'excel'),
            scroller = TRUE,
            scrollX = TRUE,
            scrollY = "700px",
            scrollResize = TRUE,
            deferRender = TRUE,
            autoWidth = TRUE
          )
        )
        
      })
        

      output$biplot <- renderPlot({
        res <- get_pcaX()
        ph <- colorby()
        shiny::req(res, ph)
        par(cex = 1.2, mar = c(6, 4, 4, 0.8))
        pcaexplorer.plot_biplot(
          res, res$Y, xpc = 1, ypc = 2,
          colorby_var = ph, arrows = arrowtype(),
          arrow.col = "black", arrow.lwd = 3,
          numarrows = n_arrows(),
          var.scale = vscale(),
          cex = 1.2, cex.labels = 1,
          legend = TRUE,
          add.ellipse = ellipse()
        )
      })

      output$biplot_2d <- plotly::renderPlotly({
        res <- get_pcaX()
        ph <- colorby()
        shiny::req(res, ph)
        pcaexplorer.plot_biplot_2d(
          res, res$Y, xpc = 1, ypc = 2,
          colorby_var = ph, arrows = arrowtype(),
          arrow.col = "black",
          numarrows = n_arrows(),
          var.scale = vscale(),
          legend = TRUE,
          add.ellipse = ellipse()
        )
      })

      output$biplot_3d <- plotly::renderPlotly({
        res <- get_pcaX()
        ph <- colorby()
        shiny::req(res, ph)
        pcaexplorer.plot_biplot_3d(
          res, res$Y, xpc = 1, ypc = 2, zpc = 3,
          colorby_var = ph, arrows = arrowtype(),
          arrow.col = "black",
          numarrows = n_arrows(),
          var.scale = vscale(),
          legend = TRUE,
          add.ellipse = ellipse()
        )
      })

      output$variance_proportion <- renderPlot({
        res <- get_pcaX()
        shiny::req(res)
        pcaexplorer.plot_variance_proportion(res)
      })

      output$variance_cumulative <- renderPlot({
        res <- get_pcaX()
        shiny::req(res)
        pcaexplorer.plot_variance_cumulative(res)
      })

    }
  )
}
