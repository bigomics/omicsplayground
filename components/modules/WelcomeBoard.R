##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

WelcomeBoard <- function(id, auth, enable_upload, r_global) {
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

    observeEvent(input$btn_example_data, {
      r_global$load_example_trigger <- r_global$load_example_trigger + 1
    })

    observeEvent(input$btn_upload_data, {
      if(enable_upload) {
        bigdash.openSidebar()
        bigdash.selectTab( session, "upload-tab" )
      } else {
        shinyWidgets::sendSweetAlert(
          session = session,
          title = "Upload disabled",          
          text ='Sorry, upload of new data is disabled for this account.',
          type = "warning",
          btn_labels = "OK",
          closeOnClickOutside = FALSE
        )
      }
    })

    observeEvent(input$btn_load_data, {
      # close the right sidebar
      #shinyjs::runjs("$('#settings-container').trigger('click');")
      #shinyjs::runjs("$('#settings-container').trigger('mouseleave');")

      bigdash.openSettings(lock=TRUE)
      shinyjs::runjs("$('#settings-container').trigger('mouseenter');")      
      bigdash.openSidebar()

      bigdash.selectTab( session, "load-tab" )
    })

  })
}

WelcomeBoardInputs <- function(id) {
  return(NULL)
}

WelcomeBoardUI <- function(id) {
  ns <- shiny::NS(id) ## namespace
  div(
    id = "welcome-page",
    div(
      class = "row",
      style = "min-height:540px;height:60vh;",
      id = "welcome-content",      
      div(
        class = "col-md-12",
        br(),
        br(),
        div(shiny::textOutput(ns("welcome")), id = "welcome-text"),
        h2("What would you like to do today?"),
        br(),
        br(),
        br()
      )
    ),
    div(
      class = "row",
      style = "max-height:35vh;padding:20px 0 20px 0;vertical-align:bottom;",
      id = "welcome-buttons",
      div(
        class = "col-md-5",
        h3("I am new..."),
        shiny::actionButton(
          ns("btn_example_data"),
          label = "Load example dataset",
          class = "btn btn-outline-info welcome-btn"
        )
      ),
      div(
        class = "col-md-7",
        h3("I'm an existing user..."),
        shiny::actionButton(
          ns("btn_upload_data"),
          label = "Upload new data",
          class = "btn btn-outline-info welcome-btn"
        ),
        shiny::actionButton(
          ns("btn_load_data"),
          label = "Use my saved data",
          class = "btn btn-outline-primary welcome-btn"
        )
      )
    ),
    br()
  )
}
