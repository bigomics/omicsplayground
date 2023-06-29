#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
     pgx_rl = reactiveVal()
     pgx <- reactive({
        browser()
        req(input$pgx_path)
        load(normalizePath(input$pgx_path))
        pgx_rl <- pgx

     })

    server <- TcgaBoard('tcga', pgx_rl)
}
