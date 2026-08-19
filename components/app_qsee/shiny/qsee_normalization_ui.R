
qsee_normalization_inputs <- function(id) {
  ns <- shiny::NS(id)
  
  bigdash::tabSettings(
    shiny::selectInput(ns("colorby"), "Color by:", choices = NULL),
    shiny::sliderInput(
      ns("noise"),
      "Noise amount:",
      min = 0,
      max = 2,
      value = 1,
      step = 0.1
    )
  )
}

qsee_normalization_ui <- function(id) {
  ns <- shiny::NS(id)
  OmicsBoardUI(
    id = ns("board"),
    title = "Normalization",
    bigdash::bd_visibility_probe(ns),
    shiny::uiOutput(ns("ui_output"), class = "html-fill-item html-fill-container")
  )
}

qsee_normalization_ui_output <- function(ns) {
  box_infotext <- "Boxplots of log-expression values (random feature subset) grouped by phenotype. Ideal normalization removes technical variation so that phenotype groups align across methods."
  hist_infotext <- "Density distributions of log-intensities per sample after each normalization method. Good normalization reduces batch effects while preserving biological signal."
  pca_infotext <- "PCA of samples after each normalization (colored by selected phenotype). Effective normalization should cluster samples by biology rather than technical factors."

  bslib::navset_tab(
    bslib::nav_panel(
      title = "Boxplots",
      bslib::layout_columns(
        col_widths = 12,
        height = "calc(100vh - 140px)",
        row_heights = list("auto", 1),
        bs_alert("Compare different normalization methods on boxplots, densities and PCA. Includes noise injection for testing robustness."),
        bslib::layout_columns(
          height = "100%",
          PlotModuleUI(
            ns("box_plots"),
            title = "Boxplots by normalization method",
            info.text = box_infotext,
            caption = box_infotext,
            plotlib = "plotly",
            height = c("100%", "80vh")
          )
        )
      )
    ),
    bslib::nav_panel(
      title = "Histograms",
      bslib::layout_columns(
        col_widths = 12,
        height = "calc(100vh - 140px)",
        PlotModuleUI(
          ns("histograms"),
          title = "Density histograms",
          info.text = hist_infotext,
          caption = hist_infotext,
          plotlib = "plotly",
          height = c("100%", "80vh")
        )
      )
    ),
    bslib::nav_panel(
      title = "PCA",
      bslib::layout_columns(
        col_widths = 12,
        height = "calc(100vh - 140px)",
        PlotModuleUI(
          ns("pca_plots"),
          title = "PCA by normalization method",
          info.text = pca_infotext,
          caption = pca_infotext,
          plotlib = "plotly",
          options = shiny::tagList(
            shiny::checkboxInput(ns("show_labels"), "Show sample labels", value = FALSE)
          ),
          height = c("100%", "80vh")
        )
      )
    )
  )
}

