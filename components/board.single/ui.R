#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {

    board = options()$board
    
    print(board)
    print(board)
    print(board)
    board = "Tcga"
    source('../app/R/global.R')

    root_opg <- get_opg_root()

    print(root_opg)
    
    source(file.path(root_opg,'components/golem_utils/app_config.R'))
    source(file.path(root_opg,'components/golem_utils/run_app.R'))
    source(file.path(root_opg,'components/golem_utils/run_dev.R'))

    directory <- file.path(root_opg, glue::glue('components/board.{board}/R/'))  # Specify the directory path
    file_paths <- list.files(directory, full.names = TRUE)  # Get the full file paths in the directory
    lapply(file_paths, source)  # Source all the files in the directory

    ui_files <- list_files_safe(path = file.path(root_opg,'components/ui/'))

    for (ui_file in ui_files) {
        source(file.path(root_opg,'components/ui/', ui_file))
    }

    # attempt to evaluate board functions
    
    # board_inputs <- grep("inputs", ls(), value = TRUE, ignore.case = TRUE)

    # header
    
    header <- shiny::tagList(
        tags$head(
        golem::favicon(),
        golem::bundle_resources(
        path = file.path(root_opg,'components/app/R/www'),
        app_title = "dev_board"
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
                board,
                paste0(board,"-tab")
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
                paste0(board,"-tab"),
                TcgaInputs("tcga"), # eval(parse(text = paste0(board_inputs, glue::glue('("{board}")'))))
                TcgaUI("tcga")
            )
        )
    )

}