#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  
  # attemp 1

      trigger_server <- reactive({
        req(input$pgx_path)
        load(input$pgx_path)
        server <- TcgaBoard('tcga', pgx)
      })
      observeEvent(input$pgx_path, {
          trigger_server()
      })


  # attemp 2

    # loaded_pgx <- reactive({
    #   req(input$pgx_path)
    #   load(input$pgx_path)
    #   names(pgx)
    #   return(pgx)
    # })
  
    # server <- TcgaBoard('tcga', loaded_pgx())

  # working version

      # load("../../../data/example-data.pgx")
      # server <- TcgaBoard('tcga', pgx)

  
}