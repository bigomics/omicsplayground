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
  boards <- boards[!boards %in% c("pathway","connectivity","enrichment","featuremap","intersection", "pcsf", "signature","wgcna")]

  message(boards)

  # create folder snap is not available
  if(!dir.exists("snap")){
    dir.create("snap")
  }


  AppDriverLog <- lapply(boards, function(board){
    # get error from AppDriver and save it as error_log
    message(board)
    #board = "dataview"
    #board = boards[1]
    try(AppDriver$stop(), silent = TRUE)
    AppDriver = NULL
    AppDriver <- tryCatch(
      {
        shinytest2::AppDriver$new(
          normalizePath("../../components/dev"),
          timeout = 10000,
          options = list(
            board = board,
            use_example_data = FALSE,
            shiny.error = function(e) {
                return("Error in shiny.error")

              }
          ),
          shiny_args = list(port = 8080)
        )},
        error = function(e) {
              # append error log to error_logs list 
              return(e)
            }
      )

    if(class(AppDriver)[[1]] == "rlang_error"){
      AppDriver$message
      return(AppDriver$message)
    }

    loadPGX <- tryCatch({
      AppDriver  # Display the AppDriver object
      # basic info for the app
      AppDriver$get_url()
      AppDriver$get_logs()

      # get pgx file path
      pgx_file <- normalizePath("../../data/example-data.pgx")

      pgx_file

      # get input/output values
      AppDriver$get_value(input = "pgx_path")
      AppDriver$get_values()

      # update pgx_path to an actual path
      AppDriver$set_inputs("pgx_path" = pgx_file)

      # check that the path is updated
      AppDriver$get_value(input = "pgx_path")
      # get screenshot to the folder snap
            
      AppDriver$get_screenshot(file.path("snap", paste0(board, ".png")))
      
      # use the get_logs to check if we have any error
      df <- data.frame(AppDriver$get_logs())
      AppDriver$stop()
      return(df)
    },
    error = function(e) {
      error = NULL
      try(error <- data.frame(AppDriver$get_logs()), silent = TRUE)
      
      if(is.null(error)){
        error <- conditionMessage(e)
      }
      # return message from error
      return(error)
    })
    AppDriver$stop()

    return(loadPGX)

  })

  # check if any board has error, via error log

  boards_with_error <- boards[sapply(AppDriverLog, function(x) any(grepl("Error in", x$message)))]

  # pass if length of boards_with_error is 0
  expect_equal(length(boards_with_error), 0)
})