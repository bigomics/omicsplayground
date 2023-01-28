#' ##
#' ## This file is part of the Omics Playground project.
#' ## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
#' ##
#'
#' #' UI code for table code: expression board
#' #'
#' #' @param id
#' #' @param label
#' #' @param height
#' #' @param width
#' #'
#' #' @export
#' compare_table_corr_score_ui <- function(id,
#'                                         label='',
#'                                         height,
#'                                         width) {
#'
#'   ns <- shiny::NS(id)
#'
#'   tableWidget(ns("table"))
#'
#' }
#'
#' #' Server side table code: expression board
#' #'
#' #' @param id
#' #' @param watermark
#' #'
#' #' @export
#' compare_table_corr_score_server <- function(id,
#'                                             watermark=FALSE){
#'   moduleServer( id, function(input, output, session) {
#'
#'     score_table.RENDER <- shiny::reactive({
#'
#'       #code here
#'
#'     })
#'
#'     score_table_info = ""
#'
#'     score_table <- shiny::callModule(
#'       tableModule, id = "table",
#'       func = "", ## ns=ns,
#'       info.text = table_info,
#'       title = tags$div(
#'         HTML('')),
#'       height = c(,),
#'       width = c(,)
#'     )
#'     return(score_table)
#'   }
#'   )
#' }