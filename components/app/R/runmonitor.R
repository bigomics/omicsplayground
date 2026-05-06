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
    col_widths = c(-2, 7, -3),
    style = "height: min(90%,700)",
    fill = TRUE,
    shiny::div(id = "navheader-current-section", HTML("Runs")),
    p("Monitor and inspect the details of computation runs"),
    bslib::card(DT::DTOutput(ns("runtable")))
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
      df <- data.frame(
        job = c("first","second"),
        start = rep(as.character(Sys.time()),2),
        duration = c("0:35 hours","-"),
        status = c("completed","running")
      )
      DT::datatable(df)
    })
    
  }) 
}
