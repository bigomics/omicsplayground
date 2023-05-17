#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  # Your application server logic
  load('data/example-data.pgx')

  diffexpr <- ExpressionBoard('diffexpr',
                              pgx = pgx)
  enrich <- EnrichmentBoard('enrich',
                            pgx = pgx,
                            selected_gxmethods = diffexpr$selected_gxmethods)
  server <- IntersectionBoard('isect',
                              pgx = pgx,
                              selected_gxmethods = diffexpr$selected_gxmethods,
                              selected_gsetmethods = enrich$selected_gsetmethods)
}
