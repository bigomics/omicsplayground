# test single board minimal components

### board specific files ###

error_logs <- list()

# get error from AppDriver and save it as error_log
options = list(board = "tcga")

AppDriver <- tryCatch(
  {
     shinytest2::AppDriver$new(
      normalizePath("components/board.single"),
      options = list(board = "tcga")
      )
  },
  error = function(e) {
    # append error log to error_logs list
    print(e)
  }
)

AppDriver  # Display the AppDriver object
# basic info for the app
AppDriver$get_url()
AppDriver$get_logs()

# get pgx file path
pgx_file <- normalizePath("data/example-data.pgx")

pgx_file

# get input/output values
AppDriver$get_value(input = "pgx_path")
AppDriver$get_values()


# update pgx_path to an actual path
AppDriver$set_inputs("pgx_path" = pgx_file)

# check that the path is updated
AppDriver$get_value(input = "pgx_path")
AppDriver$getEventLog()
# get error log

# get stderror for the app
AppDriver$get_logs()

# exit the app
AppDriver$finalize()
AppDriver$snapshot()

shiny::runApp("components/board.single/")

