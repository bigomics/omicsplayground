#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request, resources, path) {
    
    setwd(path)

    source(glue::glue('components/board.tcga/dev/app_ui.R'))
    
    
    # header
    
    header <- shiny::tagList(
        tags$head(
        golem::favicon(),
        golem::bundle_resources(
        path = file.path(path,'components/app/R/www'),
        app_title = "test"
        ),
        shinyjs::useShinyjs(),
        sever::useSever(),
        bigLoaders::addBigLoaderDeps(),
        firebase::useFirebase(firestore = TRUE, analytics = TRUE)
    )
    )

    # sidebar
    sidebar <- bigdash::sidebar(
        "Menu",
        bigdash::sidebarMenu(
            "Signature",
            bigdash::sidebarMenuItem(
                "TCGA survival (beta)",
                "tcga-tab"
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
        textInput("pgx_path", NULL),
        settings = bigdash::settings(
            "Settings"
        ),
        bigdash::bigTabs(
            bigdash::bigTabItem(
                "tcga-tab",
                TcgaInputs("tcga"),
                TcgaUI("tcga")
            )
        )
    )

}