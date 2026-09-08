setwd('components/app/R')
source('global.R')
source('ui.R')
source('server.R')
options(shiny.otel.collect = "none") 

# Accept optional port argument (default 3838). Usage: Rscript dev/run_app_headless.R [PORT]
args <- commandArgs(trailingOnly = TRUE)
port <- if (length(args) > 0) as.integer(args[1]) else 3838

shiny::runApp(
    appDir = ".",
    launch.browser = FALSE,
    host = "0.0.0.0",
    port = port
)
