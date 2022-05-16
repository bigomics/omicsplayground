##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

WelcomeBoard <- function(id, auth) {}
WelcomeBoardInputs <- function(id) {}
WelcomeBoardUI <- function(id) {}


WelcomeBoard <- function(id, auth)
{
  moduleServer(id, function(input, output, session) 
  {
    ns <- session$ns ## NAMESPACE

    output$welcome <- shiny::renderText({
        name <- auth$name()
        dbg("[HomeBoard] name = ",name)        
        if(name %in% c("",NA,NULL)) {
          welcome <- "Welcome back..."
        } else {
          first.name <- strsplit("ivo kwee",split="[@ .]")[[1]][1]
          first.name <- paste0(toupper(substring(first.name,1,1)),
                               substring(first.name,2,nchar(first.name)))
          welcome <- paste0("Welcome back ",first.name,"...")
        }
        welcome
    })

    observeEvent( input$load_data, {
      ## shinyjs::click("")
    })

    observeEvent( input$upload_new, {
      ## shinyjs::click("")
    })
       
  })
}

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
        tags$button(
          id = "init-load-data",
          "Use pre-loaded data",
          class = "btn btn-outline-primary welcome-btn"
        ),
        "Click to load a previously uploaded dataset.", placement="bottom"),
      withTooltip(
        tags$a(
          id = "init-upload-data",
          "Upload new data",
          class = "btn btn-outline-primary welcome-btn"
        ),
        "Click to upload some new data",
        placement="bottom"
      )
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
