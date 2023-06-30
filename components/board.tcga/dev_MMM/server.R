#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  # attempt 1, using reactive container
  
  loaded_pgx <- reactive({
    req(input$pgx_path)
    load(input$pgx_path)
    names(pgx)
    return(pgx)
    }
  )
  
  # loading from absolute path works, we want to pass a reactval to modify this path
  
  # load("C:/code//omicsplayground/data/example-data.pgx")

  server <- TcgaBoard('tcga', loaded_pgx())
}