
PCAexplorer_ui <- function(id) {
  ns <- shiny::NS(id)

  mod_info <- shiny::HTML("<h4>PCA explorer</h4> Explore PCA components.\n")

  ## Static sidebar widgets (never destroyed → no flicker / no sticky state).
  ## Board-dependent visibility uses a CSS class toggled from the server
  ## (toggleClass only — never jQuery show/hide, which re-inits ionSliders
  ## and flashes the panel). Nested arrow opts use client-side
  ## conditionalPanel on the checkbox so that path never hits the server.
  ##
  ## Per-board:
  ##   Biplot:        colorby, arrowtype, show_arrows, numarrows, var_scale, ellipse
  ##   PCA pairs:     colorby, ellipse
  ##   Loadings:      (none)
  ##   Feature PCA:   show_arrows, numarrows, var_scale, min_fc
  bslib::page_fillable(
    title = "PCA explorer",
    padding = 0,
    ## safe if main app already called useShinyjs(); needed for test app
    shinyjs::useShinyjs(),
    shiny::tags$style(shiny::HTML(sprintf(
      paste0(
        "#%s .pca-sb-hide { display: none !important; }",
        "#%s .bslib-sidebar-layout > .sidebar { overflow-y: auto; }"
      ),
      ns("pca_root"), ns("pca_root")
    ))),
    shiny::div(
      id = ns("pca_root"),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 220,
          position = "right",
          open = "always",
          mod_info,
          shiny::div(
            id = ns("sb_colorby"),
            shiny::selectInput(ns("colorby"), "Color by:", choices = NULL)
          ),
          shiny::div(
            id = ns("sb_show_arrows"),
            style = "margin-bottom: -20px;",
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
              ),
              hr()
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
        ),
        bslib::navset_tab(
          id = ns("board"),
          bslib::nav_panel(
            title = "Biplot",
            value = "Biplot",
            PCAexplorer_biplot_ui(ns("biplot"))
          ),
          bslib::nav_panel(
            title = "PCA pairs",
            value = "PCA pairs",
            PCAexplorer_pairs_ui(ns("pairs"))
          ),
          bslib::nav_panel(
            title = "Loadings",
            value = "Loadings",
            PCAexplorer_loadings_ui(ns("loadings"))
          ),
          bslib::nav_panel(
            title = "Feature PCA",
            value = "Feature PCA",
            PCAexplorer_feature_ui(ns("feature"))
          )
        )
      )
    )
  )
}
