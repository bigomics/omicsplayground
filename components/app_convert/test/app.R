##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Minimal standalone Shiny app to test the Gene ID Converter applet
## (app_convert) in isolation, outside the full Playground application.
##
## Run with:
##   Rscript -e 'shiny::runApp("components/app_convert/test")'
## or open this file in RStudio and click "Run App".
##

library(shiny)
library(bslib)
library(DT)

source("../R/convert_ui.R", encoding = "UTF-8")
source("../R/convert_server.R", encoding = "UTF-8")

ui <- convert_ui("convert")

server <- function(input, output, session) {
  convert_server("convert")
}

shinyApp(ui, server)
