test_that("example data loads with no error",{
  # source aux functions
  source("aux-test-functions.R")
  
  # test single board minimal components

  # get all board names
  boards <- list.dirs(path = "../../components", full.names = FALSE, recursive = FALSE)
  # split "." and get second name
  boards <- sapply(strsplit(boards, split = "\\."), function(x) x[2])
  boards <- boards[!is.na(boards)]

  # remove upload, loading and user from boards
  boards <- boards[!boards %in% c("upload", "loading", "user")]

  # remove problematic boards
  boards <- boards[!boards %in% c("tcga", "signature")]

  authentication <- options()$authentication

  AppLog <- lapply(boards, function(board) {
    # get error from App and save it as error_log
    message(board)
    # board = "wordcloud"
    # board = boards[1]
    App <- shinytest2::AppDriver$new(
      normalizePath("../../dev/board.launch"),
      timeout = 35000,
      height = 1080,
      width = 1920,
      seed = 2910,
      variant = shinytest2::platform_variant(),
      options = list(
        board = board,
        authentication = authentication,
        use_example_data = FALSE
      ),
      shiny_args = list(port = 8080)
    )
    App$get_values(input = TRUE)

    withr::defer(App$stop())

    pgx_file <- normalizePath("../../data/mini-example/example-data-mini.pgx")
    App$set_inputs("pgx_path" = pgx_file)
    if(board == "enrichment") {
      App$set_inputs("enrichment-gs_fdr" = 0.5)
      App$wait_for_idle(duration = 10000, timeout = 50000)
    }
    tabs <- searchTabs(board)
    if (!is.null(tabs)){
      lapply(tabs, function(tab){
        App$run_js(generate_js_click_code(tab))
        if(board == "connectivity") {
          duration <- 1000000
          App$wait_for_idle(duration = 10000, timeout = duration)
        } else if (board == "clustering") {
          duration <- 50000
          App$wait_for_idle(duration = 15000, timeout = duration)
        } else {
          duration <- 50000
          App$wait_for_idle(duration = 3000, timeout = duration)
        }
        
        App$expect_screenshot(cran = TRUE, name = paste0(board, "_", tab), threshold = 10, selector = "viewport")
      })
    } else {
      App$wait_for_idle(duration = 3000)
      App$expect_screenshot(cran = TRUE, name = board, threshold = 10, selector = "viewport")
    }
  })
})
