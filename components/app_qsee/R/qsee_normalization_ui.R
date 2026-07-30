

qsee_normalization_ui <- function(id) {
  ns <- shiny::NS(id)
  
  mod_info <- HTML("<h4>Normalization analysis</h4> Compare different normalization methods.\n")

  normalize_panel <- bslib::page_sidebar(
    height = "calc(100vh - 100px)",
    sidebar = bslib::sidebar(
      position = "right",
      open = "always",
      mod_info,
      selectInput(ns("colorby"), "Color by:", choices = NULL),
      shiny::sliderInput(
        ns("noise"),
        "Noise amount:",
        min = 0,
        max = 2,
        value = 1,
        step = 0.1
      )
    ),
    bslib::layout_columns(
      col_widths = 12,
      ##row_heights = c(1,1),
      height="100%",
      bslib::navset_tab(
        bslib::nav_panel(
          title = "Boxplots",
          plotOutput(ns("box_plots"), height="700px")
        ),
        bslib::nav_panel(
          title = "Histograms",
          plotOutput(ns("histograms"), height="700px")
        ),
        bslib::nav_panel(
          title = "PCA",
          plotOutput(ns("pca_plots"), height="700px")
        )
      )
    )
  )
}

    
