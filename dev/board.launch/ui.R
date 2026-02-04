##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

app_ui <- function() {
    board = options()$board
    use_example_data = options()$use_example_data

    # handle pgx file input
    pgx_file <- ""

    pgx_file_opt <- options()$pgx_file
    if (!is.null(pgx_file_opt) && nchar(pgx_file_opt) > 0) {
        pgx_file <- pgx_file_opt
    } else if (!is.null(use_example_data) && use_example_data == TRUE) {
        pgx_file <- file.path(OPG, "data/example-data.pgx")
    }

    # handle board inputs and dependencies

    directory <- file.path(OPG, glue::glue('components/board.{board}/R/'))  # Specify the directory path
    file_paths <- list.files(directory, full.names = TRUE, pattern="*[.][rR]$")  

    for (file_path in file_paths) {
        source(file_path)
    }

    ui_files <- list.files(path = file.path(OPG,'components/ui/'), pattern="*[.][rR]$")
    for (ui_file in ui_files) {
        source(file.path(OPG,'components/ui/', ui_file))
    }

    board_input <- grep("inputs", ls(envir = .GlobalEnv), value = TRUE, ignore.case = TRUE)
    # find the function name
    board_input <- board_input[grepl(board, board_input, ignore.case = TRUE)]

    # check if the length matches the function
    length <- nchar(board) + nchar("inputs")
    
    board_input <- board_input[which(length == nchar(board_input))]

    board_input_fn <- get(board_input)

    # to the same for ui
    
    ui_fn_name <- glue::glue("{board}ui")
    board_ui <- grep(ui_fn_name, ls(envir = .GlobalEnv), value = TRUE, ignore.case = TRUE)
    length <- nchar(board) + nchar("ui")
    
    board_ui <- board_ui[which(length == nchar(board_ui))]

    board_ui_fn <- get(board_ui)

    # header

     header <- shiny::tagList(
      #shiny::tags$head(htmltools::includeHTML("www/hubspot-embed.js")),
      ##    gtag2, ## Google Tag Manager???
      shiny::tags$head(shiny::tags$script(src = "static/temp.js")),
      shiny::tags$head(shiny::tags$script(src = "static/test_trigger.js")),
      shiny::tags$head(shiny::tags$script(src = "static/dropdown-extra.js")),
      shiny::tags$head(shiny::tags$script(src = "static/copy-info-helper.js")),
      shiny::tags$head(shiny::tags$script(src = "static/add-tick-helper.js")),
      shiny::tags$head(shiny::tags$link(rel = "stylesheet", href = "static/styles.min.css")),
      shiny::tags$head(shiny::tags$link(rel = "shortcut icon", href = "static/favicon.ico")),
      shinyjs::useShinyjs(),
      waiter::use_waiter(),
      sever::useSever(),
      bigLoaders::addBigLoaderDeps(),
      #firebase::useFirebase(firestore = TRUE, analytics = TRUE),
      #shinybusy::busy_start_up(
      #  text = tags$h2("\nPrepping your personal playground..."), mode = "auto",
      #  background = "#2780e3", color = "#ffffff",
      #  loader = shinybusy::spin_epic("hollow-dots", color = "#FFF")
      #)
    )

    # sidebar
    sidebar <- bigdash::sidebar(
        "Menu",
        bigdash::sidebarMenu(
            "Single Board",
            bigdash::sidebarMenuItem(
                board,
                paste0(board,"-tab")
            ),
            bigdash::sidebarMenuItem(
                "Upload PGX",
                paste0("pgx-tab")
                )
        )
    )

   big_theme2 <- bigdash::big_theme()
    big_theme2 <- bslib::bs_add_variables(big_theme2,
      "grid-breakpoints" = "map-merge($grid-breakpoints, ('xxxl': 2400px))",
      .where = "declarations"
    )

    bigdash::bigPage(
        header,
        shiny.i18n::usei18n(i18n),
        title = "Omics Playground v3",
        theme = big_theme2,
        sidebar = sidebar,
        navbar = bigdash::navbar(
            tags$img(
                id = "logo-bigomics",
                src = "static/bigomics-logo.png",
                width = "110"
            ),
            center = tags$div(
                shiny::div(shiny::textOutput("current_dataset"), class = "current-dataset"),
                shiny::div(shiny::textOutput("current_user"), class = "current-user")
        )

            
        ),
        settings = bigdash::settings(
            "Settings"
        ),
        bigdash::bigTabs(
            bigdash::bigTabItem(
                paste0(board,"-tab"),
                board_input_fn(board),
                board_ui_fn(board)
            ),
        bigdash::bigTabs(
            bigdash::bigTabItem(
                paste0("pgx-tab"),
                board_input_fn("hello"),
                textInput("pgx_path", label = NULL, value = pgx_file, placeholder = "Absolute path to pgx object", width = "100%")
            )
            
        )
        )
    )
}
