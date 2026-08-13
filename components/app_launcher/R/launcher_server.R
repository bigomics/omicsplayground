##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' The application server-side logic
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @export
launcher_server <- function(id, parent) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    shiny::observeEvent(input$logo_click, {
      ui.showAboutModal()
    })

    observeEvent(input$launch_playground, {
      bslib::nav_select("app-sidebar", "Dashboard", session=parent)
    })

    observeEvent(input$launch_across, {
      bslib::nav_select("app-sidebar", "AcrossDatasets", session=parent)
    })

    observeEvent(input$launch_prism, {
      bslib::nav_select("app-sidebar", "Prism", session=parent)
    })
    
    observeEvent(input$launch_qsee, {
      bslib::nav_select("app-sidebar", "Qsee", session=parent)
    })

    observeEvent(input$launch_methylome, {
      bslib::nav_select("app-sidebar", "Methylome", session=parent)
    })

    observeEvent(input$launch_idconvert, {
      bslib::nav_select("app-sidebar", "IDconvert", session=parent)
    })

  })
}
