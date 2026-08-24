##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Minimal standalone Shiny app to test the IDAT Converter applet
## (app_idat) in isolation, outside the full Playground application.
##
## Run with:
##   Rscript -e 'shiny::runApp("components/app_idat/test")'
## or open this file in RStudio and click "Run App".
##

library(shiny)
library(bslib)
library(DT)
library(future)
library(promises)

## IDAT files run to tens of MB per channel; Shiny's 5MB default rejects them.
## 999 MB matches components/app/R/global.R, so the standalone app rejects the
## same uploads the real one does.
options(shiny.maxRequestSize = 999 * 1024^2)

## read_idats() takes minutes; the full app sets this in global.R:141.
future::plan(future::multisession)

source("../R/idat_ui.R", encoding = "UTF-8")
source("../R/idat_server.R", encoding = "UTF-8")

ui <- idat_ui("idat")

server <- function(input, output, session) {
  idat_server("idat")
}

shinyApp(ui, server)
