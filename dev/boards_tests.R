# test single board minimal components

### board specific files ###

# get error from driver and save it as error_log
error_log <- list()

driver <- tryCatch(
      {
        shinytest::ShinyDriver$new(
          path = "components/board.tcga/dev_MMM/"
        )
      },
      error = function(e) {
        error_log = e
      }
    )

# basic info for the app
driver$getUrl()
driver$getDebugLog()
driver$listWidgets()$input
driver$listWidgets()$output

# get pgx file path
pgx_file <- normalizePath("data/example-data.pgx")

pgx_file

# I added a textInput for pgx_path
driver$getValue("pgx_path")

# update pgx_path to an actual path
driver$setValue("pgx_path", pgx_file)

# check that the path is updated
driver$getValue("pgx_path")

# get error log

# get stderror for the app
driver$getDebugLog()
driver$listWidgets()$output
driver$getValue("error_log")

# driver is still running even if we stop the app
driver
# exit the app
driver$finalize()


shiny::runApp("components/board.tcga/dev_MMM/")