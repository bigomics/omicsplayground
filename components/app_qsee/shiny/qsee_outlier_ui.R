

qsee_outlier_ui <- function(id) {
  ns <- shiny::NS(id)
  
  mod_info <- HTML("<h4>Missing value analysis</h4> Analyze imputation methods and missing value distribution.\n")

  normalize_panel <- bslib::page_sidebar(
    sidebar = bslib::sidebar(
      position = "right",
      open = "always",
      mod_info,
      shiny::selectInput(ns("colorby"), "Color by:", choices=NULL)
    ),
    bslib::navset_tab(
      bslib::nav_panel(
        title = "Z-scores",
        bslib::card(
          plotly::plotlyOutput(ns("outlier_zscores"), height="720px")
        )
      ),
      bslib::nav_panel(
        title = "Heatmap",
        bslib::layout_columns(
          col_widths = c(6,6),
          bslib::card(iheatmapr::iheatmaprOutput(ns("outlier_heatmap"), height="720px")),
          bslib::card(plotly::plotlyOutput(ns("outlier_pca"), height="720px"))
        )
      )
    )
  )

  shiny::tagList(
    qsee_visibility_probe(ns),
    normalize_panel
  )
}

