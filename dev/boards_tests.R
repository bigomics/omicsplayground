# test single board minimal components

library(golem)
library(playbase)
options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode
options(shiny.port = httpuv::randomPort())
golem::detach_all_attached()

ui_files <- list_files_safe(path = 'components/ui/')

for (ui_file in ui_files) {
  source(file.path('components/ui/', ui_file))
}

### board specific files ###

driver <- shinytest::ShinyDriver$new(
  path = "components/board.tcga/dev_MMM/"
  )

driver$getUrl()
driver$getDebugLog()
driver$listWidgets()$input
driver$listWidgets()$output

  # loadTimeout = NULL,
  # checkNames = TRUE,
  # debug = "shiny_console",
  # phantomTimeout = 50000,
  # seed = NULL,
  # cleanLogs = FALSE,
  # shinyOptions = list(),
  # renderArgs = NULL,
  # options = list()
  # )


pgx_file <- normalizePath("data/example-data.pgx")

driver$setValue("pgx_path", pgx_file)

driver$getValue("tcga-contrast")

# get all values for the tcga-contrast input
driver$getAllValues("tcga-contrast")

driver$finalize()