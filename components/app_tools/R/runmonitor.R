##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


RunMonitorUI <- function(id) {
  ns <- shiny::NS(id)
  
  ui <- bslib::layout_columns(
    col_widths = c(-1, 10, -1),
    style = "height: min(90%,700)",
    fill = TRUE,
    shiny::div(id = "navheader-current-section", HTML("Runs")),
    h4("Monitor and inspect the details of computation runs"),
    bslib::layout_columns(
      col_widths = c(9, 3),
      bslib::card(
        bslib::card_header("Submitted runs"),
        bslib::card_body(DT::DTOutput(ns("runtable")))
      ),
      bslib::card(
        bslib::card_header("Run info"),
        bslib::card_body(textOutput(ns("info")))
      )
    )
  )

  return(ui)
}

#' The application server-side logic
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @export
#' Note: pgx needs to be reactiveValues
#'
#'
RunMonitorServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## ----------------- plot output -------------------------------
    output$runtable <- DT::renderDT({
      tt <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      df <- data.frame(
        job = c("my_first_run","my_second_run"),
        start = rep(tt,2),
        duration = c("0:35 hours","-"),
        status = c("completed","running")
      )
      DT::datatable(df)
    })

    output$info <- renderText({
      paste("Information")
    })
    
  }) 
}
