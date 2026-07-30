## Feature PCA board: feature-space PCA colored by trait correlation,
## PC–trait correlation heatmap, and compact loading barplots.

PCAexplorer_feature_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_columns(
    col_widths = c(7, 5),
    navset_card_tab(
      full_screen = TRUE,
      title = "Feature PCA",
      bslib::nav_spacer(),
      nav_panel("plotly", plotly::plotlyOutput(ns("feature_pca"), height = "800px")),
      nav_panel("3D", plotly::plotlyOutput(ns("feature_pca_3d"), height = "800px")),
      nav_panel("all", shiny::plotOutput(ns("feature_pca_all"), height = "800px"))
    ),
    bslib::layout_columns(
      col_widths = 12,
      row_heights = c(0.7, 1),
      bslib::card(full_screen = TRUE, shiny::plotOutput(ns("trait_correlation"), height = "400px")),
      bslib::card(full_screen = TRUE, shiny::plotOutput(ns("loadings2"), height = "400px"))
    )
  )
}

PCAexplorer_feature_module <- function(id, get_pcaX, get_normX, rY,
                                       numarrows = shiny::reactive(10),
                                       show_arrows = shiny::reactive(TRUE),
                                       min_fc = shiny::reactive(0),
                                       var_scale = shiny::reactive(0.65)) {
  shiny::moduleServer(
    id,
    function(input, output, session) {

      n_arrows <- shiny::reactive({
        if (isTRUE(show_arrows())) {
          n <- suppressWarnings(as.integer(numarrows()))
          if (!is.finite(n) || n < 1L) 10L else n
        } else {
          0L
        }
      })

      fc_cut <- shiny::reactive({
        v <- suppressWarnings(as.numeric(min_fc()))
        if (!is.finite(v) || v < 0) 0 else v
      })

      vscale <- shiny::reactive({
        v <- suppressWarnings(as.numeric(var_scale()))
        if (!is.finite(v) || v <= 0) 0.65 else v
      })

      output$feature_pca <- plotly::renderPlotly({
        res <- get_pcaX()
        normX <- get_normX()
        shiny::req(res, normX)
        pcaexplorer.plot_feature_pca(
          res, normX = normX, Y = rY(),
          numarrows = n_arrows(),
          min_fc = fc_cut(),
          var.scale = vscale()
        )
      })

      output$feature_pca_3d <- plotly::renderPlotly({
        res <- get_pcaX()
        normX <- get_normX()
        shiny::req(res, normX)
        pcaexplorer.plot_feature_pca_3d(
          res, normX = normX, Y = rY(),
          numarrows = n_arrows(),
          min_fc = fc_cut(),
          var.scale = vscale()
        )
      })

      output$feature_pca_all <- shiny::renderPlot({
        res <- get_pcaX()
        normX <- get_normX()
        shiny::req(res, normX)
        pcaexplorer.plot_feature_pca_all(
          res, normX = normX, Y = rY(),
          min_fc = fc_cut(),
          var.scale = vscale(),
          show_arrows = isTRUE(show_arrows())
        )
      })

      output$trait_correlation <- shiny::renderPlot({
        res <- get_pcaX()
        shiny::req(res)
        pcaexplorer.plot_trait_correlation(res, min_fc = fc_cut())
      })

      output$loadings2 <- shiny::renderPlot({
        res <- get_pcaX()
        shiny::req(res)
        pcaexplorer.plot_loadings_bars_compact(
          res, npc = 2, ntop = 40,
          normX = get_normX(), Y = rY(),
          min_fc = fc_cut()
        )
      })

    }
  )
}
