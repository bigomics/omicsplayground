#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  options(shiny.error = function() {
      # The error message is on the parent environment, it is
      # not passed to the function called on error
      parent_env <- parent.frame()
      error <- parent_env$e
      error_log(error)
    })

  
  # attempt 1
  error_log <- reactiveVal("No error")
  
  trigger_server <- reactive({
    req(input$pgx_path)
    load(input$pgx_path)
    server <- TcgaBoard('tcga', pgx)
  })
  
  observeEvent(input$pgx_path, {
    trigger_server()
  })
  
  # Render error log as output text
  output$error_log <- renderText({
    logs <- error_log()
    if (is.list(logs)) {
      logs <- paste(logs, collapse = "\n")
    }
    pgx_logs
    }
  )
  
}