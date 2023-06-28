# test single board minimal components

library(golem)
library(playbase)
options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode
options(shiny.port = httpuv::randomPort())
golem::detach_all_attached()
source('components/app/R/global.R')
source('components/golem_utils/app_config.R')
source("components\\golem_utils\\run_dev.R")
source("components\\golem_utils\\run_app.R")

board = "board.tcga"
load('data/example-data.pgx')


ui_files <- list_files_safe(path = 'components/ui/')

for (ui_file in ui_files) {
  source(file.path('components/ui/', ui_file))
}
source('components/golem_utils/app_config.R')
source('components/golem_utils/run_app.R')
source('components/golem_utils/run_dev.R')

### board specific files ###

source(glue::glue('components/{board}/dev/app_ui.R'))
# source(glue::glue('components/{board}/dev/app_server.R'))


app_server <- function(input, output, session) {

    server <- TcgaBoard('tcga', pgx)
}


r_files <- list_files_safe(path = glue::glue('components/{board}/R'))

for (r_file in r_files) {
  source(file.path(glue::glue('components/{board}/R/'),r_file))
}

# Run the application
run_app(options=options)
