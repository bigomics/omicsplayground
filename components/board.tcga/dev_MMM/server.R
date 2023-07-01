#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  
  # attempt 1
  error_log <- reactiveVal("No error")
  
  trigger_server <- reactive({
    req(input$pgx_path)
    load(input$pgx_path)
    server <- TcgaBoard('tcga', pgx)
  })
  
  observeEvent(input$pgx_path, {
    tryCatch(
      {
        trigger_server()
      },
      error = function(e) {
        error_log(e)
      }
    )
  })
  
  # Render error log as output text
  output$error_log <- renderText({
  logs <- error_log()
  if (is.list(logs)) {
    logs <- paste(logs, collapse = "\n")
  }
  logs
})
  
}