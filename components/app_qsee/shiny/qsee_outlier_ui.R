

qsee_outlier_inputs <- function(id) {
  ns <- shiny::NS(id)

  bigdash::tabSettings(
    shiny::selectInput(ns("colorby"), "Color by:", choices=NULL)
  )
}

qsee_outlier_ui <- function(id) {
  ns <- shiny::NS(id)

  zscore_infotext <- "Z-score based outlier detection. Multi-panel barplots show per-sample and per-feature outlier scores (or summary statistics). Horizontal lines mark common significance thresholds."
  heatmap_infotext <- "Clustered heatmap of the most outlier-ish features (row-centered and scaled). Dendrograms highlight extreme patterns."
  pca_infotext <- "PCA plot of samples colored by selected phenotype. Outliers often appear as distant points."

  outlier_panel <-
    bslib::navset_tab(
      bslib::nav_panel(
        title = "Z-scores",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 140px)",
          PlotModuleUI(
            ns("outlier_zscores"),
            title = "Z-score outliers",
            info.text = zscore_infotext,
            caption = zscore_infotext,
            plotlib = "plotly",
            height = c("100%", "80vh")
          )
        )
      ),
      bslib::nav_panel(
        title = "Heatmap & PCA",
        bslib::layout_columns(
          col_widths = c(6, 6),
          height = "calc(100vh - 140px)",
          bslib::card(
            full_screen = TRUE,
            iheatmapr::iheatmaprOutput(ns("outlier_heatmap"), height = "calc(100vh - 260px)")
          ),
          PlotModuleUI(
            ns("outlier_pca"),
            title = "PCA with outliers",
            info.text = pca_infotext,
            caption = pca_infotext,
            plotlib = "plotly",
            options = shiny::tagList(
              shiny::checkboxInput(ns("show_labels"), "Show sample labels", value = TRUE)
            ),
            height = c("100%", "80vh")
          )
        )
      )
    )

  OmicsBoardUI(
    id = ns("board"),
    title = "Outlier analysis",
    qsee_visibility_probe(ns),
    outlier_panel
  )
}

