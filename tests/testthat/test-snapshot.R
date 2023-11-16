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
  boards <- boards[!boards %in% c("tcga")]

  AppLog <- lapply(boards, function(board){
    # get error from App and save it as error_log
    message(board)
    #board = "wordcloud"
    #board = boards[1]
    App <- shinytest2::AppDriver$new(
      normalizePath("../../dev/board.launch"),
      timeout = 20000,
      height = 1080,
      width = 1920,
      seed = 2910,
      variant = shinytest2::platform_variant(),
      options = list(
        board = board,
        use_example_data = FALSE
      ),
      shiny_args = list(port = 8080)
    )
    
    withr::defer(App$stop())
    
    pgx_file <- normalizePath("../../data/mini-example/example-data-mini.pgx")
    App$set_inputs("pgx_path" = pgx_file)

    App$wait_for_idle(duration=10000)
    
    # App$expect_values(cran = TRUE) # TODO: file bug about this...
    App$expect_screenshot(cran = TRUE, name = board, threshold = 10, selector = "viewport")
  })
})
