
library(optparse)

# Define the options
option_list <- list(
  make_option(c("-d", "--data"), type="character", help="Path to pgx data files")
)

# Parse the command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$data)) {
  stop("Please provide a path to data files using the -d or --data option.")
}

#opt <- list()
#opt$data <- "/mnt/c/code/omicsplayground/data/mini-example"
pgx_files <- list.files(path = opt$data,  full.names = TRUE, recursive = TRUE)

# get all files that end in pgx
pgx_files <- unique(pgx_files[grepl("pgx$", pgx_files)])

# add absolute path
pgx_files <- normalizePath(pgx_files)

# get pgx file name
pgx_file_name <- basename(pgx_files)

# test single board minimal components

# get all board names
boards <- list.dirs(path = "components", full.names = FALSE, recursive = FALSE)
# split "." and get second name
boards <- sapply(strsplit(boards, split = "\\."), function(x) x[2])
boards <- boards[!is.na(boards)]

# remove upload, loading and user from boards
boards <- boards[!boards %in% c("upload", "loading", "user")]

# remove problematic boards
boards <- boards[!boards %in% c("featuremap","intersection", "wordcloud")]


# df with name as pgx_name and each column will result of be a board

results <- list()

# test single board minimal components
for (pgx_file in pgx_files) {
  #pgx_file <- pgx_files[1]
  message(pgx_file)

  AppLog <- lapply(boards, function(board){
    # get error from App and save it as error_log
    message(board)
    #board = "dataview"
    #board = boards[7]
    App <- tryCatch(
      {
        shinytest2::AppDriver$new(
          normalizePath("dev/board.launch"),
          timeout = 20000,
          options = list(
            board = board,
            shiny.error = function(e) {
                return("Error in shiny.error")

              }
          ),
          shiny_args = list(
            port = 8080
            )
        )
        },
        error = function(e) {
              # append error log to error_logs list 
              return(e)
            }
      )
    withr::defer(App$stop())
    if(class(App)[[1]] == "rlang_error"){
      App$message
      return(App$message)
    }

    loadPGX <- tryCatch({
      App  # Display the App object
      # basic info for the app
      App$get_url()
      App$get_logs()

      # get input/output values
      App$get_value(input = "pgx_path")
      App$get_values()

      # update pgx_path to an actual path
      App$set_inputs("pgx_path" = pgx_file)

      # check that the path is updated
      App$get_value(input = "pgx_path")
      # get error log

      # use the get_logs to check if we have any error
      df <- data.frame(App$get_logs())
      App$stop()
      return(df)
    },
    error = function(e) {
      error = NULL
      try(error <- data.frame(App$get_logs()), silent = TRUE)
      
      if(is.null(error)){
        error <- conditionMessage(e)
      }
      # return message from error
      return(error)
    })
    return(loadPGX)
})

# check if any board has error, via error log

boards_with_error <- sapply(AppLog, function(x) any(grepl("Error in", x$message)))

results[[pgx_file]] <- boards_with_error
}

pgx_check_results <- data.frame(do.call(rbind, results))

colnames(pgx_check_results) <- boards

# save results to opt$data folder
write.csv(pgx_check_results, file = paste0(opt$data, "pgx_check_results.csv"))
