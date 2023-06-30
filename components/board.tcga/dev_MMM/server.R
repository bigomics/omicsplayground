#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  source('C:/code/omicsplayground/components/app/R/global.R')

  source('C:/code//omicsplayground/components/golem_utils/app_config.R')
  source('C:/code//omicsplayground/components/golem_utils/run_app.R')
  source('C:/code//omicsplayground/components/golem_utils/run_dev.R')
  board = "board.tcga"
  ui_files <- list_files_safe(path = 'C:/code//omicsplayground/components/ui/')

  for (ui_file in ui_files) {
      source(file.path('C:/code//omicsplayground/components/ui/', ui_file))
  }

  r_files <- list_files_safe(path = normalizePath('C:/code//omicsplayground/components/board.tcga/R'))

  for (r_file in r_files) {
      source(file.path(glue::glue('C:/code/omicsplayground/components/{board}/R/'),r_file))
  }

  load("C:/code//omicsplayground/data/example-data.pgx") # this somehow does not work with

  server <- TcgaBoard('tcga', pgx)
}