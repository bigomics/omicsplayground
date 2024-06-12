##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

DatasetReportUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::actionButton(
    ns("action"), "Generate Report",
    width = "auto", class = "quick-button"
  )
}

DatasetReportServer <- function(
    id,
    auth,
    callbackR = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    ## email text input validator
    
    showModal <- function() {
      body <- tagList(
        HTML("<center><h3><b>and earn some swag!</b></h3><p><p>Invite your friends to Omics Playground and earn some exclusive Bigomics swag like cool stickers, our 'Friendly Monster' T-Shirt or one of our awesome sustainable Dopper water bottles. Read more about it <a href='https://bigomics.ch/invite' target='_blank'><u>here</u></a>.<br><br>"),
        div(
          shiny::textInput(
            inputId = ns("available_datasets"),
            label = "", placeholder = "Email address..."
          ),
          style = "margin-top: -30px;"
        ),
        shiny::actionButton(ns("generate_report_action"), "Submit", class = "btn btn-primary")
      )

      modal <- shiny::modalDialog(
        title = NULL,
        bsutils::modalHeader(
          div(class = "modal-title", "Create a report"),
          style = "background-color: #f0f9fd;"
        ),
        body,
        footer = NULL,
        size = "l",
        easyClose = TRUE,
        tags$style(".modal-dialog {width: 720px;}"),
        tags$style(".modal-content {background-color: #f0f9fd;}"),
        tags$style(".modal-header {padding: 0px;}")
      )

      shiny::showModal(modal)
    }

    r_click <- shiny::reactiveVal(0)
    ext_click <- function() {
      r_click(r_click() + 1)
    }

    click <- shiny::reactive({
      r_click() + input$action
    })

    shiny::observeEvent(
      {
        list(r_click(), input$action)
      },
      {
        if (r_click() || input$action) {
          showModal()
        }
      }
    )

    shiny::observeEvent(input$generate_report_action, {
      
    })

    list(
      click = ext_click ## exported function!
    )
  }) ## end of moduleServer
}
