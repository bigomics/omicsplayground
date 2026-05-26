##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' The application server-side logic
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @export
tools_server <- function(id, parent) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$runtool1, {
      bslib::nav_select("app-sidebar", "Prism", session=parent)
    })

    
  })
}
