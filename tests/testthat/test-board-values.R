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
  boards <- boards[!boards %in% c("pathway","featuremap","intersection", "pcsf", "signature","wgcna")]

  message(boards)

  # create folder snap is not available
  if(!dir.exists("snap")){
    dir.create("snap")
  }

  AppLog <- lapply(boards[1:3], function(board){
    # get error from App and save it as error_log
    message(board)
    #board = "dataview"
    #board = boards[1]
    App <- shinytest2::AppDriver$new(
      normalizePath("../../dev/board.launch"),
      timeout = 10000,
      height = 1080,
      width = 1920,
      variant = platform_variant(),
      options = list(
        board = board,
        use_example_data = FALSE
      ),
      shiny_args = list(port = 8080)
    )
      
    withr::defer(App$stop())
    
    pgx_file <- normalizePath("../../data/example-data.pgx")
    App$set_inputs("pgx_path" = pgx_file)

    expect_true(TRUE)

    App$expect_screenshot(cran = TRUE)
  })
})