#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # Your application server logic
  pgx <- playbase::pgx.load("data/example-data.pgx")
  server <- PcsfBoard("pcsf", pgx)
}
