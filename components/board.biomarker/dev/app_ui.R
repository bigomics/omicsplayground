#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {

    # header
    header <- shiny::tagList(
        golem_add_external_resources('board.biomarker')
    )

    # sidebar
    sidebar <- bigdash::sidebar(
        "Menu",
        bigdash::sidebarMenu(
            "Signature",
            bigdash::sidebarItem(
                "Find biomarkers",
                "bio-tab"
            )
        )
    )

    big_theme <- bslib::bs_add_variables(
        bigdash::big_theme(),
        "grid-breakpoints" = "map-merge($grid-breakpoints, ('xxxl': 2400px))",
        .where = "declarations"
    )

    bigdash::bigPage(
        header,
        title = "Omics Playground v3",
        theme = big_theme,
        sidebar = sidebar,
        navbar = bigdash::navbar(
            tags$img(
                id = "logo-bigomics",
                src = "assets/img/bigomics.png",
                width = "110",
            )
        ),
        settings = bigdash::settings(
            "Settings"
        ),
        bigdash::bigTabs(
            bigdash::bigTabItem(
                "bio-tab",
                BiomarkerInputs("bio"),
                BiomarkerUI("bio")
            )
        )
    )

}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(board) {
  golem::add_resource_path(
    "www",
    'components/app/R/www'
    #app_sys("app/www")
  )

  tags$head(
    golem::favicon(),
    golem::bundle_resources(
      path = 'components/app/R/www',
      #path = app_sys("app/www"),
      app_title = board
    ),
    shinyjs::useShinyjs(),
    sever::useSever(),
    bigLoaders::addBigLoaderDeps(),
    firebase::useFirebase(firestore = TRUE, analytics = TRUE)
  )
}
