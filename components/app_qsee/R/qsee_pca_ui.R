

qsee_pca_ui <- function(id) {
  ns <- shiny::NS(id)
  
  mod_info <- HTML("<h4>PCA explorer</h4> Explore PCA components.\n")

  normalize_panel <- bslib::page_fillable(
    padding = 0,
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width = 220,
        position = "right",
        open = "always",
        mod_info,
        selectInput(ns("arrows"),"Arrows:",choices=c("loadings","pheno")),
        selectInput(ns("maxarrows"),"Max. arrows:", choices = c(3,5,10,20),
          selected=10),
      ),
      bslib::navset_tab(
        bslib::nav_panel(
          title = "PCA pairs",
          bslib::layout_columns(
            col_widths = c(7,5),
            bslib::card(plotOutput(ns("pairs"),height="800px")),
            bslib::layout_columns(
              col_widths = 12,
              row_heights = c(1,1),
              bslib::card(plotOutput(ns("variance"), height="400px")),
              bslib::card(plotOutput(ns("trait_correlation"), height="400px"))
            )
          )
        ),
        bslib::nav_panel(
          title = "Biplot",
          bslib::layout_columns(
            col_widths = c(6,6),
            bslib::card(plotOutput(ns("biplot"),height="800px")),
            bslib::card(plotOutput(ns("biplot2"),height="800px"))
          )
        ),      
        bslib::nav_panel(
          title = "Feature PCA",
          bslib::layout_columns(
            col_widths = c(7,5),
            bslib::card(plotOutput(ns("feature_pca"),height="800px")),
            br()
          )
        )      
      )
    )
  )
}

    
