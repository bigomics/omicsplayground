# test single board minimal components

library(golem)
library(playbase)
options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode
options(shiny.port = httpuv::randomPort())
golem::detach_all_attached()
source('components/app/R/global.R')

board = "board.tcga"
load('data/example-data.pgx')

source('components/golem_utils/app_config.R')
source('components/golem_utils/run_app.R')
source('components/golem_utils/run_dev.R')


# root dir


ui_files <- list_files_safe(path = 'components/ui/')

for (ui_file in ui_files) {
  source(file.path('components/ui/', ui_file))
}

### board specific files ###

source(glue::glue('components/board.tcga/dev_MMM/app_ui.R'))

app_server <- function(input, output, session) {



    source('C:/Users/Xavier/OneDrive/BigOmics/GitHub/omicsplayground/components/app/R/global.R')
    source('C:/Users/Xavier/OneDrive/BigOmics/GitHub/omicsplayground/components/golem_utils/app_config.R')
    source('C:/Users/Xavier/OneDrive/BigOmics/GitHub/omicsplayground/components/golem_utils/run_app.R')
    source('C:/Users/Xavier/OneDrive/BigOmics/GitHub/omicsplayground/components/golem_utils/run_dev.R')
    board = "board.tcga"
    ui_files <- list_files_safe(path = 'C:/Users/Xavier/OneDrive/BigOmics/GitHub/omicsplayground/components/ui/')

    for (ui_file in ui_files) {
        source(file.path('C:/Users/Xavier/OneDrive/BigOmics/GitHub/omicsplayground/components/ui/', ui_file))
    }

    r_files <- list_files_safe(path = normalizePath('C:/Users/Xavier/OneDrive/BigOmics/GitHub/omicsplayground/components/board.tcga/R'))

    for (r_file in r_files) {
        source(file.path(glue::glue('C:/Users/Xavier/OneDrive/BigOmics/GitHub/omicsplayground/components/{board}/R/'),r_file))
    }

  load("C:/Users/Xavier/OneDrive/BigOmics/GitHub/omicsplayground/data/example-data.pgx") # this somehow does not work with

  server <- TcgaBoard('tcga', pgx)
}

r_files <- list_files_safe(path = normalizePath('C:/Users/Xavier/OneDrive/BigOmics/GitHub/omicsplayground/components/board.tcga/R'))

for (r_file in r_files) {
  source(file.path(glue::glue('components/{board}/R/'),r_file))
}

onStart = NULL
enableBookmarking = NULL
uiPattern = "/"
resources <- golem_add_external_resources("board.tcga")

app = shinyApp(
  ui = app_ui(resources = resources, path = getwd()),
  server = app_server,
  onStart = onStart,
  options = options,
  enableBookmarking = enableBookmarking,
  uiPattern = uiPattern
)

driver <- shinytest::ShinyDriver$new(
  path = app
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
  driver$listWidgets()$output

  driver$getAllValues("TcgaBoard")

  pgx_file <- normalizePath("data/example-data.pgx")

  driver$setValue("pgx_path", pgx_file)

  driver$getValue("pgx_path")

  driver$finalize()
