

qsee_imputation_inputs <- function(id) {
  ns <- shiny::NS(id)

  bigdash::tabSettings(
    shiny::selectInput(ns("colorby"), "Color by:", choices=NULL),
    div(
      shiny::sliderInput(ns("marlevel"), "Add MAR:", 0, 1, 0, 0.05),
      shiny::sliderInput(ns("mnarlevel"), "Add MNAR:", 0, 1, 0, 0.05)
    )
  )
}

qsee_imputation_ui <- function(id) {
  ns <- shiny::NS(id)
  OmicsBoardUI(
    id = ns("board"),
    title = "Missing value analysis",
    bigdash::bd_visibility_probe(ns),
    shiny::uiOutput(ns("ui_output"), class = "html-fill-item html-fill-container")
  )
}

qsee_imputation_ui_output <- function(ns) {
  dist_infotext <- "Missingness patterns versus intensity: barplots of missing ratio per intensity bin, feature-wise scatterplots, histograms of missing ratios, retention curves, and per-sample missing rates."
  hist_infotext <- "Density distributions of observed versus imputed intensity values for each imputation method."
  pca_infotext <- "PCA of samples after each imputation method (samples colored by phenotype)."
  valid_infotext <- "Imputation validation: imputed values versus known (simulated MAR/MNAR) true values. Higher Pearson correlation indicates better imputation performance."

  bslib::navset_tab(
      bslib::nav_panel(
        title = "Distribution",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 140px)",
          PlotModuleUI(
            ns("distributions"),
            title = "Missingness distributions",
            info.text = dist_infotext,
            caption = dist_infotext,
            plotlib = "plotly",
            height = c("100%", "80vh")
          )
        )
      ),
      bslib::nav_panel(
        title = "Histogram",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 140px)",
          PlotModuleUI(
            ns("histograms"),
            title = "Observed vs imputed densities",
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
            title = "PCA after imputation",
            info.text = pca_infotext,
            caption = pca_infotext,
            plotlib = "plotly",
            options = shiny::tagList(
              shiny::checkboxInput(ns("show_labels"), "Show sample labels", value = FALSE)
            ),
            height = c("100%", "80vh")
          )
        )
      ),
      bslib::nav_panel(
        title = "Heatmap",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 140px)",
          bslib::card(
            full_screen = TRUE,
            bslib::card_header("Imputed data correlation heatmaps"),
            qsee_plotly_hm_grid_ui(ns, "heatmap", n = 6L, ncol = 3L, height = "calc(50vh - 130px)")
          )
        )
      ),
      bslib::nav_panel(
        title = "Validation",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 140px)",
          PlotModuleUI(
            ns("validation_scatter"),
            title = "Imputation accuracy (simulated missing)",
            info.text = valid_infotext,
            caption = valid_infotext,
            plotlib = "plotly",
            height = c("100%", "80vh")
          )
        )
      )
    )
}

