# test single board minimal components

### board specific files ###

driver <- shinytest::ShinyDriver$new(
  path = "components/board.tcga/dev_MMM/"
  )


# when I run the same thing over runApp, and place C:\\code\\omicsplayground\\data\\example-data.pgx in text input, the pgx gets loaded ;)
shiny::runApp("components/board.tcga/dev_MMM/")

# basic info for the app
driver$getUrl()
driver$getDebugLog()
driver$listWidgets()$input
driver$listWidgets()$output

# get pgx file path
pgx_file <- normalizePath("data/example-data.pgx")

# I added a textInput for pgx_path
driver$getValue("pgx_path")

# update pgx_path to an actual path
driver$setValue("pgx_path", pgx_file)

# check that the path is updated
driver$getValue("pgx_path")

# exit the app
driver$finalize()
