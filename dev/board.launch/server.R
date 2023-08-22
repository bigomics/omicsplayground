#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  # list functions in global
  
  board = options()$board

  server_fn_name <- glue::glue("{board}board")
  board_server <- grep(server_fn_name, ls(envir = .GlobalEnv), value = TRUE, ignore.case = TRUE)
  length <- nchar(board) + nchar("board")
    
  board_server <- board_server[which(length == nchar(board_server))]

  board_server_fn <- get(board_server)

  trigger_server <- reactive({
        req(input$pgx_path)
        browser()
        pgx <- playbase::pgx.load(input$pgx_path)
        server <- board_server_fn(board, pgx = pgx)

  })
  
  observeEvent(input$pgx_path, {
    trigger_server()
  })
}