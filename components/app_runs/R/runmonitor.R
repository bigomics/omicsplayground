##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


RunMonitorUI <- function(id) {
  ns <- shiny::NS(id)
  
  ui <- bslib::layout_columns(
    col_widths = c(-2, 8, -2),
    row_heights = list("auto","auto",1),
    style = "height: min(90%,700); margin-bottom: 10px;",
    fill = TRUE,
    shiny::div(id = "navheader-current-section", HTML("Runs")),
    h4("Monitor and inspect the details of computation runs"),
    bslib::layout_columns(
      col_widths = c(8, 4),
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
      DT::datatable(
        df,
        plugins = "scrollResize",
        fillContainer = TRUE,
        options = list(
          dom = "t",
          pageLength = 999,
          scrollResize = TRUE
        ),
        class = "compact hover"
      )
    })

    output$info <- renderText({
      paste("Information")
    })
    
  }) 
}
