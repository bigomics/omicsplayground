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


ui_files <- list_files_safe(path = 'components/ui/')

for (ui_file in ui_files) {
  source(file.path('components/ui/', ui_file))
}

### board specific files ###

source(glue::glue('components/board.tcga/dev_MMM/app_ui.R'))

app_server <- function(input, output, session) {
  r_files <- list.files("C:\\code\\omicsplayground\\components\\board.tcga\\R",full.names = TRUE, include.dirs = TRUE)

  for (r_file in r_files) {
    source(r_file)
  }
  
  pgx_rl <- reactiveVal(NULL)

  observeEvent(input$pgx_path, {
    req(input$pgx_path)
    load(normalizePath(input$pgx_path))
    
    load(normalizePath(driver$getValue("pgx_path")))
    
    message("pgx loaded")
    
    pgx_rl(pgx)
  })

  server <- TcgaBoard('tcga', pgx_rl())
}

r_files <- list_files_safe(path = normalizePath('components/board.tcga/R')

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
  path = app,
  loadTimeout = NULL,
  checkNames = TRUE,
  debug = "shiny_console",
  phantomTimeout = 50000,
  seed = NULL,
  cleanLogs = FALSE,
  shinyOptions = list(),
  renderArgs = NULL,
  options = list()
  )

  driver$listWidgets()$input
  driver$listWidgets()$output

  driver$getAllValues("TcgaBoard")

  pgx_file <- normalizePath("data/example-data.pgx")

  driver$setValue("pgx_path", pgx_file)

  driver$getValue("pgx_path")

  driver$finalize()