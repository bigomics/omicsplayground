##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

WelcomeBoardInputs <- function(id) {
  ## ns <- shiny::NS(id)  ## namespace
  ## bigdash::tabSettings(
  ##   shiny::actionLink(ns("module_info"), "Tutorial", icon = shiny::icon("youtube"))
  ## )
  return(NULL)
}

WelcomeBoardUI <- function(id) {
  ns <- shiny::NS(id)  ## namespace

  fillCol(
    id = "welcome-page",
    style = "text-align:center;background-color:#eaf7fd;",
    flex = c(1,NA,NA,NA,1,2),
    height = "95vh",

    br(),
    div(shiny::textOutput(ns("welcome")), id="welcome-text"),
    h3("What would you like to do today?", style="padding: 10px;"),
    
    fillRow(
      id="welcome-buttons",
      height = 150,
      withTooltip(
        shiny::actionButton(
          ns("load_data"), label="Use pre-loaded data",
          class="btn btn-outline-primary welcome-btn"),
        "Click to load a previously uploaded dataset.", placement="bottom"),
      withTooltip( actionButton(
        ns("upload_new"), label="Upload new data",
        class="btn btn-outline-primary welcome-btn")
       ,"Click to upload some new data", placement="bottom")
    ),
    
    div(
      id="welcome-subtext",
      HTML("<B>BigOmics Playground. Never Stop Discovering.</B><br>
          BigOmics is focused on one thing — helping life scientists see and understand their omics
          data. Our mission is to create smart tools and make advanced omics analysis accessible to
          everyone. Want to know more? Read our paper \"Omics Playground: a comprehensive self-
          service platform for visualization, analytics and exploration of Big Omics Data\".")),
    br()
    
  ) ## end of fill-col
}
