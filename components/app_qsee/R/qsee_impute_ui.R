

qsee_imputation_ui <- function(id) {
  ns <- shiny::NS(id)
  
  mod_info <- HTML("<h4>Missing value analysis</h4> Analyze imputation methods and missing value distribution.\n")

  normalize_panel <- bslib::page_sidebar(
    sidebar = bslib::sidebar(
      position = "right",
      open = "always",
      mod_info,
      shiny::selectInput(ns("colorby"), "Color by:", choices=NULL),
      div(
        shiny::sliderInput(ns("marlevel"), "Add MAR:", 0, 1, 0, 0.05),
        shiny::sliderInput(ns("mnarlevel"), "Add MNAR:", 0, 1, 0, 0.05)
      )
    ),
    bslib::navset_tab(
      bslib::nav_panel(
        title = "Distribution",
        bslib::card(
          plotOutput(ns("distributions"), height="720px")
        )
      ),
      bslib::nav_panel(
        title = "Histogram",
        bslib::card(
          plotOutput(ns("histograms"), height="720px")
        )
      ),
      bslib::nav_panel(
        title = "PCA",
        bslib::card(
          plotOutput(ns("pca_plots"), height="720px")
        )
      ),
      bslib::nav_panel(
        title = "Heatmap",
        bslib::card(
          plotOutput(ns("heatmap_plots"), height="720px")
        )
      ),
      bslib::nav_panel(
        title = "Validation",
        bslib::card(
          plotOutput(ns("validation_scatter"), height="720px")
        )
      )
    )
  )
}

    
