##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics Sagl. All rights reserved.
##

WelcomeBoard <- function(id, auth) {}
WelcomeBoardInputs <- function(id) {}
WelcomeBoardUI <- function(id) {}


WelcomeBoard <- function(id, auth, r_global) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    output$welcome <- shiny::renderText({
      name <- auth$name()
      dbg("[WelcomeBoard] name =",name)
      if (name %in% c("", NA, NULL)) {
        welcome <- "Welcome back..."
      } else {
        first.name <- strsplit(name, split = "[@ .]")[[1]][1]
        first.name <- paste0(
          toupper(substring(first.name, 1, 1)),
          substring(first.name, 2, nchar(first.name))
        )
        welcome <- paste0("Welcome back ", first.name, "...")
      }
      welcome
    })

    observeEvent(input$init_example_data, {
      r_global$load_example_trigger <- TRUE
    })

    observeEvent(input$init_upload_data, {
      bigdash.openSidebar()
      bigdash.selectTab( session, "upload-tab" )
    })

    observeEvent(input$init_load_data, {
      # close the right sidebar
      shinyjs::runjs("$('#settings-container').trigger('click');")
      shinyjs::runjs("$('#settings-container').trigger('mouseleave');")

      bigdash.openSidebar()
      bigdash.selectTab( session, "load-tab" )
    })

  })
}

WelcomeBoardInputs <- function(id) {
  return(NULL)
}

WelcomeBoardUI <- function(id, enable_upload=TRUE) {
  ns <- shiny::NS(id) ## namespace

  ## if upload enabled show button, otherwise empty
  upload_button <- "  "
  if(enable_upload) {
    upload_button <- shiny::actionButton(
      ns("init_upload_data"),
      label = "Upload new data",
      class = "btn btn-outline-info welcome-btn"
    )
  }

  div(
    id = "welcome-page",
    br(),
    br(),
    div(shiny::textOutput(ns("welcome")), id = "welcome-text"),
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
          ns("init_example_data"),
          label = "Try example dataset",
          class = "btn btn-outline-info welcome-btn"
        )
      ),
      div(
        class = "col-md-7",
        h3("I'm an existing user..."),
        upload_button,
        shiny::actionButton(
          ns("init_load_data"),
          label = "Use my saved data",
          class = "btn btn-outline-primary welcome-btn"
        )
      )
    ),
    br()
  )
}
