
# board_deps : vector of other board names that the main board depends on
# e.g. `launch_board('board.enrichment', c('board.expression'))` because
# board.enrichment takes input from board.expression. Any board deps need
# to also be included in the app_ui.R and app_server.R files
launch_board <- function(board, board_deps = NULL) {

  library(golem)
  library(playbase) ## install or devtools::load_all(PATH_TO_PLAYBASE)

  # Set options here
  options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode

  # Comment this if you don't want the app to be served on a random port
  options(shiny.port = httpuv::randomPort())

  # Detach all loaded packages and clean your environment
  golem::detach_all_attached()
  # rm(list=ls(all.names = TRUE))

  document_and_reload2 <-
    function (pkg = get_golem_wd(), roclets = NULL, load_code = NULL,
              clean = FALSE, export_all = FALSE, helpers = FALSE, attach_testthat = FALSE,
              ...)
    {
      #check_name_consistency(pkg)
      rlang::check_installed("pkgload")
      if (rlang::is_installed("rstudioapi") && rstudioapi::isAvailable() &&
          rstudioapi::hasFun("documentSaveAll")) {
        rstudioapi::documentSaveAll()
      }
      roxed <- try({
        golem:::roxygen2_roxygenise(package.dir = pkg, roclets = roclets,
                                    load_code = load_code, clean = clean)
      })
      if (attempt::is_try_error(roxed)) {
        golem:::cli_cat_rule("Error documenting your package")
        golem:::dialog_if_has("Alert", "Error documenting your package")
        return(invisible(FALSE))
      }
      loaded <- try({
        golem:::pkgload_load_all(pkg, export_all = export_all, helpers = helpers,
                                 attach_testthat = attach_testthat, ...)
      })
      if (attempt::is_try_error(loaded)) {
        golem:::cli_cat_rule("Error loading your package")
        golem:::dialog_if_has("Alert", "Error loading your package")
        return(invisible(FALSE))
      }
    }
  # Document and reload your package
  #document_and_reload2()

  ## SOURCE FILES ##

  # global file
  global_file <- 'components/app/R/global.R'
  source(global_file)
  #
  ## ui files from component
  ui_files <- list.files(path = 'components/ui/')
  for (ui_file in ui_files) {
    source(file.path('components/ui/', ui_file))
  }

  source('components/golem_utils/app_config.R')
  source('components/golem_utils/run_app.R')

  ### board specific files ###
  source(glue::glue('components/{board}/dev/app_ui.R'))
  source(glue::glue('components/{board}/dev/app_server.R'))
  r_files <- list.files(path = glue::glue('components/{board}/R'))
  for (r_file in r_files) {
    source(file.path(glue::glue('components/{board}/R/'),r_file))
  }

  # load any other boards
  if (!is.null(extra_boards)) {
    for (extra_board in extra_boards) {
      r_files <- list.files(path = glue::glue('components/{extra_board}/R'))
      for (r_file in r_files) {
        source(file.path(glue::glue('components/{extra_board}/R/'),r_file))
      }
    }
  }

  # Run the application
  run_app()
}


