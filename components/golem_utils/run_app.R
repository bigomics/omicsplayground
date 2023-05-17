#' Run the Shiny Application
#'
#' @param ... arguments to pass to golem_opts.
#' See `?golem::get_golem_options` for more details.
#' @inheritParams shiny::shinyApp
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(
  onStart = NULL,
  options = list(),
  enableBookmarking = NULL,
  uiPattern = "/",
  ...
) {
  golem::with_golem_options(
    app = shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = onStart,
      options = options,
      enableBookmarking = enableBookmarking,
      uiPattern = uiPattern
    ),
    golem_opts = list(...)
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(board) {
  golem::add_resource_path(
    "www",
    'components/app/R/www'
  )

  tags$head(
    golem::favicon(),
    golem::bundle_resources(
      path = 'components/app/R/www',
      app_title = board
    ),
    shinyjs::useShinyjs(),
    sever::useSever(),
    bigLoaders::addBigLoaderDeps(),
    firebase::useFirebase(firestore = TRUE, analytics = TRUE)
  )
}