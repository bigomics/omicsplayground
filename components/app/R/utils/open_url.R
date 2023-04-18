#' Open a URL using JavaScript in a Shiny App
#'
#' This function opens a given URL in a new browser window or tab using JavaScript
#' within a Shiny app. This is particularly useful when using a Docker container,
#' as the `browseURL()` function may not work as expected.
#'
#' @param url A character string representing the URL to open.
#'
#' @examples
#' \dontrun{
#' library(shiny)
#' library(shinyjs)
#'
#' ui <- fluidPage(
#'   useShinyjs(),
#'   actionButton("open_twitter", "Open Twitter")
#' )
#'
#' server <- function(input, output, session) {
#'   observeEvent(input$open_twitter, {
#'     url_to_open <- "https://twitter.com"
#'     open_url_js(url_to_open)
#'   })
#' }
#'
#' shinyApp(ui, server)
#' }
#'
#' @export

open_url_js <- function(url) {
  # Turn `'` character into %27
  url <- url%>% gsub("'", "%27", x = .)
  # Remove `\n`
  url <- url%>% gsub("\n", "", x = .)
  # Run JS code in the browser (window.open opens a new tab)
  shinyjs::runjs(
    code = paste0(
      "window.open('",
      url,
      "');"
    )
  )
}
