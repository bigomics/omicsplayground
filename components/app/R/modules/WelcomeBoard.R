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

  div(
      id = "welcome-page",
      style = "text-align:center;background-color:#eaf7fd;",
      br(),      
      br(),
      div(shiny::textOutput(ns("welcome")), id="welcome-text"),
      h3("What would you like to do today?"),
      br(),
      br(),    
      br(),
      br(),    
      div(
          class = "row",
          id = "welcome-buttons",
          div(
              class = "col-md-5",
              h4("Just want to try out?"),
              tags$a(
                  id = "init-example-data",
                  "Try example dataset",
                  class = "btn btn-outline-primary welcome-btn"
              )
          ),
          div(
              class = "col-md-7",
              h4("I'm an existing user..."),
              tags$button(
                  id = "init-load-data",
                  "Use my saved data",
                  class = "btn btn-outline-primary welcome-btn"
              ),
              tags$a(
                  id = "init-upload-data",
                  "Upload new data",
                  class = "btn btn-outline-primary welcome-btn"
              )
          )
      ),
      ## br(),
      ## div(
      ##   id="welcome-subtext",
      ##   HTML("<B>BigOmics Playground. Never Stop Discovering.</B><br>
      ##       BigOmics is focused on one thing — helping life scientists see and understand their omics
      ##       data. Our mission is to create smart tools and make advanced omics analysis accessible to
      ##       everyone. Want to know more? Read our paper \"Omics Playground: a comprehensive self-
      ##       service platform for visualization, analytics and exploration of Big Omics Data\".")),
      br()
      
  )
}
