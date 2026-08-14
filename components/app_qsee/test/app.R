##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Minimal standalone Shiny app to test the QSEE (batch-effects/QC) applet
## (app_qsee) in isolation, outside the full Playground application.
##
## The code in ../R is a work in progress and may not fully work yet.
## This harness just wires it up against example data so it can be run
## and iterated on.
##
## Run with:
##   Rscript -e 'shiny::runApp("components/app_qsee/test")'
## or ./run.sh from this directory, or open this file in RStudio and
## click "Run App".
##

library(shiny)
library(bslib)
library(dplyr)
#library(playbase)

## Make the shared styles and JS assets from the main app available.
## (matches components/app/R/global.R and ui.R)
www_path <- normalizePath("../../assets", mustWork = FALSE)
if (dir.exists(www_path)) {
  shiny::addResourcePath("custom", www_path)
  shiny::addResourcePath("assets", www_path)
} else {
  message("[test/app.R] WARNING: www dir not found at ", www_path)
}


## global `opt` list expected by the shared PlotModule UI (normally set up
## in components/app/R/global.R). Assigned into .GlobalEnv explicitly since
## shiny::runApp() evaluates this file in its own environment, and the
## sourced ui-*.R functions below look up `opt` lexically from .GlobalEnv.
assign(
  "opt",
  playbase::pgx.readOptions(file = "../../../etc/OPTIONS", default = list()),
  envir = globalenv()
)

## shared UI helpers used by qsee_ui.R / qsee_server.R (same set sourced by
## the full app in components/00SourceAll.R)
ui_files <- list.files("../../ui", pattern = "\\.R$", full.names = TRUE)
for (f in ui_files) source(f, encoding = "UTF-8")

r_files <- list.files("../R", pattern = "\\.R$", full.names = TRUE)
for (f in r_files) source(f, encoding = "UTF-8")

app_files <- list.files("../shiny", pattern = "\\.R$", full.names = TRUE)
for (f in app_files) source(f, encoding = "UTF-8")

## example dataset to drive the reactive inputs
#load("../../../data/example-data.pgx")
pgx <- playbase::pgx.load("~/Playground/omicsplayground/data/GSE10846-dlbcl-nc.pgx")
if (!exists("pgx")) pgx <- playdata::GEIGER_PGX
pgx <- NULL

ui <- bslib::page_fillable(
  ## Include the main app's styles (fonts, layout, buttons, cards, etc.)
  padding = 0, gap = 0,
  shiny::tags$head(
    shiny::tags$link(
      rel = "stylesheet",
      href = "custom/styles.min.css"
    )
  ),
  shinyjs::useShinyjs(),
  qsee_ui("qsee")
)

server <- function(input, output, session) {
  ## The qsee_server module now registers its own outputs and PlotModuleServers
  qsee_server("qsee", pgx = pgx)
}

shinyApp(ui, server)
