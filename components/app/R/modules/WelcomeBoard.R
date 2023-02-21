##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

WelcomeBoard <- function(id, auth) {}
WelcomeBoardInputs <- function(id) {}
WelcomeBoardUI <- function(id) {}


WelcomeBoard <- function(id, auth, r_global)
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

    observeEvent(input$init_example_data, {
      shinyjs::runjs("$('.tab-sidebar:eq(1)').trigger('click');")
      shinyjs::runjs("$('.sidebar-label').trigger('click');")
      r_global$load_example_trigger <- TRUE
    })

  })
}

WelcomeBoardInputs <- function(id) {
  return(NULL)
}

WelcomeBoardUI <- function(id) {
  ns <- shiny::NS(id)  ## namespace

  div(
      id = "welcome-page",
      br(),
      br(),
      div(shiny::textOutput(ns("welcome")), id="welcome-text"),
      h2("What would you like to do today?"),
      br(),
      br(),
      br(),
      div(
          class = "row",
          id = "welcome-buttons",
          div(
              class = "col-md-5",
              h3("I am new..."),
              shiny::actionButton(
                ns('init_example_data'),
                label = "Try example dataset",
                class = "btn btn-outline-info welcome-btn"
              )
          ),
          div(
              class = "col-md-7",
              h3("I'm an existing user..."),
              tags$a(
                  id = "init-upload-data",
                  "Upload new data",
                  class = "btn btn-outline-info welcome-btn"
              ),
              tags$button(
                  id = "init-load-data",
                  "Use my saved data",
                  class = "btn btn-outline-primary welcome-btn"
              )
          )
      ),
      br()

  )
}
