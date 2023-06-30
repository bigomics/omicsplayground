#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
    
    source('C:/code//omicsplayground/components/golem_utils/app_config.R')
    source('C:/code//omicsplayground/components/golem_utils/run_app.R')
    source('C:/code//omicsplayground/components/golem_utils/run_dev.R')
    source('C:/code/omicsplayground/components/app/R/global.R')
    source(file.path(getwd(),'../R/tcga_ui.R'))
    source(file.path(getwd(),'../R/tcga_server.R'))
    source(file.path(getwd(),'../R/tcga_plot_survival.R'))

    ui_files <- list_files_safe(path = 'C:/code//omicsplayground/components/ui/')

    for (ui_file in ui_files) {
        source(file.path('C:/code//omicsplayground/components/ui/', ui_file))
    }

    
    # header
    
    header <- shiny::tagList(
        tags$head(
        golem::favicon(),
        golem::bundle_resources(
        path = 'C:/code/omicsplayground/components/app/R/www',
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