##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Minimal standalone Shiny app to test the Smart Tools app-launcher UI
## (app_tools) in isolation, outside the full Playground application.
##
## Run with:
##   Rscript -e 'shiny::runApp("components/app_tools/test")'
## or ./run.sh from this directory, or open this file in RStudio and
## click "Run App".
##
## Note: clicking a tile calls bslib::nav_select("app-sidebar", ...) on
## the parent session, which only exists in the full app -- here it will
## just show a "Run" click was registered without navigating anywhere.
##
## A chooser at the top lets you instantly switch between all four UI versions:
## • Grok launcher (new rectangular mobile-style design — from scratch)
## • Claude launcher (previous app-tile grid design)
## • DS launcher (iOS/Android home-screen style with status bar)
## • Classic (original 3-card bslib layout)
##

library(shiny)
library(bslib)
library(dplyr)

source("./launcher_ui_grok.R", encoding = "UTF-8")
source("./launcher_ui_claude.R", encoding = "UTF-8")
source("./launcher_ui_classic.R", encoding = "UTF-8")
source("./launcher_ui_ds.R", encoding = "UTF-8")
source("./launcher_server.R", encoding = "UTF-8")

css_path <- "../../app/R/www/styles.min.css"
OPG = "~/Playground/omicsplayground/"
shiny::addResourcePath("www", file.path(OPG, "components/app/R/www"))

ui <- bslib::page_fluid(
  tags$head(tags$style(HTML(paste(readLines(css_path, warn = FALSE), collapse = "\n")))),
  shiny::div(
    style = "padding: 15px 20px 0 20px;",
    shiny::radioButtons(
      "ui_version", "UI version:",
      choices = c(
        "Grok launcher (new)" = "grok",
        "Claude launcher" = "claude",
        "DS launcher" = "ds",
        "Classic (pre-redesign)" = "classic"
      ),
      selected = "grok", inline = TRUE
    )
  ),
  shiny::uiOutput("chosen_ui")
)

server <- function(input, output, session) {

  output$chosen_ui <- renderUI({
    if (input$ui_version == "classic") {
      launcher_ui_classic("tools")
    } else if (input$ui_version == "claude") {
      launcher_ui_claude("tools")
    } else if (input$ui_version == "ds") {
      launcher_ui_ds("tools")
    } else if (input$ui_version == "grok") {
      launcher_ui_grok("tools")
    } else {
      # default / "grok"
      launcher_ui_grok("tools")
    }
  })

  launcher_server("tools", parent = session)
}

shinyApp(ui, server)
