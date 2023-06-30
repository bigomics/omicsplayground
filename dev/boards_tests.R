# test single board minimal components

library(golem)
library(playbase)
options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode
options(shiny.port = httpuv::randomPort())
golem::detach_all_attached()

# root dir

board = "board.tcga"


ui_files <- list_files_safe(path = 'components/ui/')

for (ui_file in ui_files) {
  source(file.path('components/ui/', ui_file))
}

### board specific files ###

source(glue::glue('components/board.tcga/dev/app_ui.R'))

app_server <- function(input, output, session) {

  source('C:/code/omicsplayground/components/app/R/global.R')

  get_opg_root <- function() {
    return("C:/code/omicsplayground")
  }
  source('C:/code//omicsplayground/components/golem_utils/app_config.R')
  source('C:/code//omicsplayground/components/golem_utils/run_app.R')
  source('C:/code//omicsplayground/components/golem_utils/run_dev.R')
  board = "board.tcga"
  ui_files <- list_files_safe(path = 'C:/code//omicsplayground/components/ui/')

  for (ui_file in ui_files) {
      source(file.path('C:/code//omicsplayground/components/ui/', ui_file))
  }

  r_files <- list_files_safe(path = normalizePath('C:/code//omicsplayground/components/board.tcga/R'))

  for (r_file in r_files) {
      source(file.path(glue::glue('C:/code/omicsplayground/components/{board}/R/'),r_file))
  }

  load("C:/code//omicsplayground/data/example-data.pgx") # this somehow does not work with

  server <- TcgaBoard('tcga', pgx)
}

r_files <- list_files_safe(path = normalizePath('C:/code//omicsplayground/components/board.tcga/R'))

for (r_file in r_files) {
  source(file.path(glue::glue('components/{board}/R/'),r_file))
}

onStart = NULL
enableBookmarking = NULL
uiPattern = "/"
resources <- golem_add_external_resources("board.tcga")


driver <- shinytest::ShinyDriver$new(
  path = "components/board.tcga/dev_MMM/"
  )

driver$getUrl()
driver$getDebugLog()
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

driver$listWidgets()$input
driver$listWidgets()$output[2]


pgx_file <- normalizePath("data/example-data.pgx")

driver$setValue("pgx_path", pgx_file)

driver$getValue("tcga-contrast")

# get all values for the tcga-contrast input
driver$getAllValues("tcga-contrast")

driver$finalize()