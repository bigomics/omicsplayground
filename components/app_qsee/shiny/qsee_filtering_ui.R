

qsee_filtering_inputs <- function(id) {
  ns <- shiny::NS(id)

  bigdash::tabSettings(
    shiny::selectInput(ns("colorby"), "Color by:", choices = NULL),
    hr(),
    shiny::sliderInput(
      ns("threshold"),
      "Threshold (top SD features):",
      min = 1,
      max = 1000,
      value = 1000,
      step = 1
    )
  )
}

qsee_filtering_ui <- function(id) {
  ns <- shiny::NS(id)
  OmicsBoardUI(
    ns,
    title = "SD Filtering",
    bigdash::bd_visibility_probe(ns),
    shiny::uiOutput(ns("ui_output"), class = "html-fill-item html-fill-container")
  )
}

qsee_filtering_ui_output <- function(ns) {
  pca_infotext <-
    "PCA of the top-N most variable features (different N per panel). Good filtering preserves major biological structure in the first few PCs."

  var_infotext <-
    "Cumulative variance explained as a function of the number of top-SD features kept. Steep initial slope shows that few high-variance features capture most signal."

  hist_infotext <-
    "Distribution of feature standard deviations across the dataset. The dashed line indicates the current threshold for top-SD features."

  bslib::navset_tab(
      bslib::nav_panel(
        title = "SD Filtering",
        bslib::layout_columns(
          col_widths = c(8, 4),
          height = "calc(100vh - 140px)",
          PlotModuleUI(
            ns("pca_vs_topsd"),
            title = "PCA (top SD features)",
            info.text = pca_infotext,
            caption = pca_infotext,
            plotlib = "plotly",
            options = shiny::tagList(
              shiny::checkboxInput(ns("show_labels"), "Show sample labels", value = FALSE)
            ),
            height = c("100%", "80vh")
          ),
          bslib::layout_columns(
            col_widths = 12,
            row_heights = c(1, 1),
            height = "calc(100vh - 140px)",
            PlotModuleUI(
              ns("variance_vs_topsd"),
              title = "Cumulative variance explained",
              info.text = var_infotext,
              caption = var_infotext,
              plotlib = "plotly",
              height = c("100%", "80vh")
            ),
            PlotModuleUI(
              ns("sd_histogram"),
              title = "SD histogram",
              info.text = hist_infotext,
              caption = hist_infotext,
              plotlib = "plotly",
              height = c("100%", "80vh")
            )
          )
        )
      )
    )
}

