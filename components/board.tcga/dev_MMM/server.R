#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  # attempt 1
  error_log <- reactiveVal(NULL)
  
  options(shiny.error = function() {
    parent_env <- parent.frame()
    error <- parent_env$e
    error_log(error)

    output$error_log <- renderText({
      logs <- error_log()
      if (is.list(logs)) {
        logs <- paste(logs, collapse = "\n")
      }
      logs
      }
    )

    return(error)
  })

  trigger_server <- reactive({
    req(input$pgx_path)
    load(input$pgx_path)
    server <- TcgaBoard('tcga', pgx)
  })
  
  observeEvent(input$pgx_path, {
    trigger_server()
  })

  # observe event error that will stop the app
  observeEvent(error_log(), {
    req(error_log())
    stop()
    shiny::stopApp()
  })
}