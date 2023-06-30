#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  load("C:/code//omicsplayground/data/example-data.pgx") # this somehow does not work with

  server <- TcgaBoard('tcga', pgx)
}