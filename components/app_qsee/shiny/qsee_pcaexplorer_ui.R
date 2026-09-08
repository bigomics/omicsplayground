##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## UI for the "PCA explorer" Qsee board (formerly the standalone
## app_pcaexplorer app, now merged in as an extra board). Thin wrapper:
## computation lives in app_qsee/R/pcaexplorer-compute.R, plotting in
## app_qsee/R/pcaexplorer-plots.R.
##

qsee_pcaexplorer_inputs <- function(id) {
  ns <- shiny::NS(id)

  ## Per-board visibility of these widgets (which sub-tab needs which
  ## widget) is toggled from the server via the "pca-sb-hide" class.
  bigdash::tabSettings(
    shiny::tags$style(shiny::HTML(".pca-sb-hide { display: none !important; }")),
    shiny::div(
      id = ns("sb_colorby"),
      shiny::selectInput(ns("colorby"), "Color by:", choices = NULL)
    ),
    shiny::div(
      id = ns("sb_show_arrows"),
      shiny::checkboxInput(ns("show_arrows"), "show arrows", value = TRUE)
    ),
    ## client-side only: no server round-trip when checkbox flips
    shiny::div(
      id = ns("sb_arrow_opts"),
      shiny::conditionalPanel(
        condition = "input.show_arrows == true",
        ns = ns,
        shiny::div(
          id = ns("sb_arrowtype"),
          shiny::selectInput(ns("arrowtype"), "Arrow type:",
            choices = c("pheno", "loadings"))
        ),
        shiny::selectInput(ns("numarrows"), "Max. arrows:",
          choices = c(3, 5, 10, 20), selected = 10),
        shiny::sliderInput(
          ns("var_scale"), "Arrow length:",
          min = 0.1, max = 1.5, value = 0.65, step = 0.05
        )
      )
    ),
    ## hidden on Biplot (default tab) until Feature PCA is selected
    shiny::div(
      id = ns("sb_min_fc"),
      class = "pca-sb-hide",
      shiny::sliderInput(
        ns("min_fc"), "Min. maxFC:",
        min = 0, max = 2, value = 0, step = 0.05
      )
    ),
    shiny::div(
      id = ns("sb_ellipse"),
      shiny::checkboxInput(ns("ellipse"), "add ellipse", value = TRUE)
    )
  )
}

qsee_pcaexplorer_ui <- function(id) {
  ns <- shiny::NS(id)
  OmicsBoardUI(
    ns,
    title = "PCA explorer",
    bigdash::bd_visibility_probe(ns),
    shiny::uiOutput(ns("ui_output"), class = "html-fill-item html-fill-container")
  )
}

qsee_pcaexplorer_ui_output <- function(ns) {
  biplot_panel <- bslib::layout_columns(
    col_widths = c(7, 5),
    height = "calc(100vh - 140px)",
    bslib::navset_card_tab(
      full_screen = TRUE,
      title = "Biplot",
      bslib::nav_spacer(),
      bslib::nav_panel("plotly", plotly::plotlyOutput(ns("biplot_2d"), height = "800px")),
      bslib::nav_panel("3D", plotly::plotlyOutput(ns("biplot_3d"), height = "800px")),
      bslib::nav_panel("table", shiny::dataTableOutput(ns("biplot_table")))
    ),
    bslib::layout_columns(
      col_widths = 12,
      row_heights = c(1, 1),
      PlotModuleUI(
        ns("variance_proportion"),
        title = "Variance vs. PC dimension",
        plotlib = "base",
        height = c("400px", "70vh")
      ),
      PlotModuleUI(
        ns("variance_cumulative"),
        title = "Variance explained",
        plotlib = "base",
        height = c("400px", "70vh")
      )
    )
  )

  pairs_panel <- bslib::layout_columns(
    col_widths = c(7, 5),
    height = "calc(100vh - 140px)",
    PlotModuleUI(
      ns("scatterpairs"),
      title = "PC scatterpair matrix",
      plotlib = "base",
      height = c("100%", "80vh")
    ),
    bslib::layout_columns(
      col_widths = 12,
      row_heights = c(1, 1),
      PlotModuleUI(
        ns("pc_boxplots"),
        title = "PC score boxplots",
        plotlib = "base",
        height = c("100%", "70vh")
      ),
      PlotModuleUI(
        ns("ftest_vs_pcdim"),
        title = "F-test vs. PC dimension",
        plotlib = "base",
        height = c("100%", "70vh")
      )
    )
  )

  loadings_panel <- bslib::layout_columns(
    col_widths = c(6, 6),
    height = "calc(100vh - 140px)",
    bslib::navset_card_tab(
      full_screen = TRUE,
      title = "Loadings",
      bslib::nav_spacer(),
      bslib::nav_panel(
        "plot",
        shiny::plotOutput(ns("loadings"), height = "800px")
      ),
      bslib::nav_panel(
        "table",
        shiny::dataTableOutput(ns("loadings_table"))
      )
    ),
    PlotModuleUI(
      ns("loadings_heatmap"),
      title = "Top-gene loading heatmap",
      plotlib = "base",
      height = c("800px", "80vh")
    )
  )

  feature_panel <- bslib::layout_columns(
    col_widths = c(7, 5),
    height = "calc(100vh - 140px)",
    bslib::navset_card_tab(
      full_screen = TRUE,
      title = "Feature PCA",
      bslib::nav_spacer(),
      bslib::nav_panel("plotly", plotly::plotlyOutput(ns("feature_pca"), height = "800px")),
      bslib::nav_panel("3D", plotly::plotlyOutput(ns("feature_pca_3d"), height = "800px")),
      bslib::nav_panel("all", shiny::plotOutput(ns("feature_pca_all"), height = "800px"))
    ),
    bslib::layout_columns(
      col_widths = 12,
      row_heights = c(0.7, 1),
      PlotModuleUI(
        ns("trait_correlation"),
        title = "PC-trait correlation",
        plotlib = "base",
        height = c("400px", "70vh")
      ),
      PlotModuleUI(
        ns("loadings2"),
        title = "Top loadings (compact)",
        plotlib = "base",
        height = c("400px", "70vh")
      )
    )
  )

  bslib::navset_tab(
    id = ns("board_tab"),
    bslib::nav_panel(title = "Biplot", value = "Biplot", biplot_panel),
    bslib::nav_panel(title = "PCA pairs", value = "PCA pairs", pairs_panel),
    bslib::nav_panel(title = "Loadings", value = "Loadings", loadings_panel),
    bslib::nav_panel(title = "Feature PCA", value = "Feature PCA", feature_panel)
  )

}
