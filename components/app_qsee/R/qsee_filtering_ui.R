

qsee_filtering_ui <- function(id) {
  ns <- shiny::NS(id)
  
  mod_info <- HTML("<h4>Missing value analysis</h4> Analyze imputation methods and missing value distribution.\n")

  normalize_panel <- bslib::page_sidebar(
    sidebar = bslib::sidebar(
      position = "right",
      open = "always",
      mod_info,
      shiny::selectInput(ns("colorby"), "Color by:", choices = NULL),
      shiny::sliderInput(
        ns("threshold"),
        "Threshold (top SD features):",
        min = 1,
        max = 1000,
        value = 1000,
        step = 1
      )
    ),
    bslib::navset_tab(
      bslib::nav_panel(
        title = "SD filtering",
        bslib::layout_columns(
          col_widths = c(7, 5),
          bslib::card(full_screen = TRUE, shiny::plotOutput(ns("pca_vs_topsd"), height = "800px")),
          bslib::layout_columns(
            col_widths = 12,
            row_heights = c(1, 1),
            bslib::card(full_screen = TRUE, shiny::plotOutput(ns("variance_vs_topsd"), height = "350px")),
            bslib::card(full_screen = TRUE, shiny::plotOutput(ns("sd_histogram"), height = "350px"))
          )
        )
      )
    )
  )
}

    
