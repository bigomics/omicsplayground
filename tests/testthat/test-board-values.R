test_that("example data loads with no error",{
  # test single board minimal components

  # get all board names
  boards <- list.dirs(path = "../../components", full.names = FALSE, recursive = FALSE)
  # split "." and get second name
  boards <- sapply(strsplit(boards, split = "\\."), function(x) x[2])
  boards <- boards[!is.na(boards)]
  
  # remove upload, loading and user from boards
  boards <- boards[!boards %in% c("upload", "loading", "user")]

  # remove problematic boards
  boards <- boards[!boards %in% c("pathway","compare","featuremap","intersection", "pcsf", "signature","wgcna")]

  AppLog <- lapply(boards, function(board){
    # get error from App and save it as error_log
    message(board)
    #board = "connectivity"
    #board = boards[1]
    App <- shinytest2::AppDriver$new(
      normalizePath("../../dev/board.launch"),
      timeout = 10000,
      height = 1080,
      width = 1920,
      variant = shinytest2::platform_variant(),
      options = list(
        board = board,
        use_example_data = FALSE
      ),
      shiny_args = list(port = 8080)
    )
      
    withr::defer(App$stop())
    
    pgx_file <- normalizePath("../../data/example-data.pgx")
    App$set_inputs("pgx_path" = pgx_file)

    # simulate a click on logo-bigomics
    App$run_js("document.getElementById('logo-bigomics').click();")

    while(App$get_value(input = "continue_test") != "TRUE"){
      Sys.sleep(100.0)
    }
    
    # App$expect_values(cran = TRUE) # TODO: file bug about this...
    App$expect_screenshot(cran = TRUE, name = board, threshold = 10)
  })
})