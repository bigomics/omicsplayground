#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  browser()
  # list functions in global
  
  board = options()$board

  server_fn_name <- glue::glue("{board}board")
  board_server <- grep(server_fn_name, ls(envir = .GlobalEnv), value = TRUE, ignore.case = TRUE)
  board_server_fn <- get(board_server)
  
  trigger_server <- reactive({
        req(input$pgx_path)
        load(input$pgx_path)
        server <- board_server_fn(board, pgx)
    
  })
  
  observeEvent(input$pgx_path, {
    trigger_server()
  })
}