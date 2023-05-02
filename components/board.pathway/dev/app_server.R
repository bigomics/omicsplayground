#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  # Your application server logic
  load('data/example-data.pgx')
  diffexpr <- ExpressionBoard('diffexpr', pgx)
  enrich <- EnrichmentBoard('enrich', pgx, diffexpr$selected_gxmethods)
  server <- FunctionalBoard('pathway', pgx, enrich$selected_gsetmethods)
}
