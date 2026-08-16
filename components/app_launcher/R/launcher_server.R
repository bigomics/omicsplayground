##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' The application server-side logic
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @export
launcher_server <- function(id, parent, load_example = NULL,
                            app_launchers = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    shiny::observeEvent(input$logo_click, {
      ui.showAboutModal()
    })

    observeEvent(input$launch_playground, {
      bslib::nav_select("app-sidebar", "Dashboard", session=parent)
    })
    
    observeEvent(input$launch_qsee, {
      run_qsee <- app_launchers[["qsee"]]
      if (!is.null(run_qsee) && is.function(run_qsee)) {
        run_qsee()
      } 
    })

    observeEvent(input$launch_across, {
      run_across <- app_launchers[["across"]]
      if (!is.null(run_across) && is.function(run_across)) {
        run_across()
      } 
      ##bslib::nav_select("app-sidebar", "AcrossDatasets", session=parent)
    })

    observeEvent(input$launch_prism, {
      bslib::nav_select("app-sidebar", "Prism", session=parent)
    })

    observeEvent(input$launch_idconvert, {
      bslib::nav_select("app-sidebar", "IDconvert", session=parent)
    })

    ## Quick action: load the example dataset and jump to the Dashboard
    observeEvent(input$load_example, {
      if (is.null(load_example)) {
        warning("[launcher_server] !!! no load_example trigger available")
        return()
      }
      if (is.null(load_example())) {
        load_example(1)
      } else {
        load_example(load_example() + 1)
      }
      bslib::nav_select("app-sidebar", "Dashboard", session = parent)
    })

    ## Quick action: go to the Upload panel
    observeEvent(input$upload_new_data, {
      bslib::nav_select("app-sidebar", "Upload", session = parent)
    })

    ## Quick action: go to the Obi AI Copilot panel
    observeEvent(input$chat_with_obi, {
      bslib::nav_select("app-sidebar", "Copilot", session = parent)
    })

  })
}
