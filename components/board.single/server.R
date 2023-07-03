#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  trigger_server <- reactive({
        req(input$pgx_path)
        load(input$pgx_path)
        server <- TcgaBoard('tcga', pgx)
    
  })
  
  observeEvent(input$pgx_path, {
    trigger_server()
  })
}