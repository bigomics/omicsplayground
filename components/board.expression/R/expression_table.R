##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' UI code for table code: expression board
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
expression_table_TableName_ui <- function(id,
                                        label='',
                                        height,
                                        width) {

  ns <- shiny::NS(id)

  tableWidget(ns("table"))

}

#' Server side table code: expression board
#'
#' @param id
#' @param watermark
#'
#' @export
expression_table_TableName_server <- function(id,
                                              watermark=FALSE){
  moduleServer( id, function(input, output, session) {

    table.RENDER <- shiny::reactive({

      #code here

    })

    table_info = ""

    score_table <- shiny::callModule(
      tableModule, id = "table",
      func = "table.RENDER", ## ns=ns,
      info.text = table_info,
      title = tags$div(
        HTML('')),
      height = c(,),
      width = c(,)
    )
    return(score_table)
  }
  )
}