#
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

UpgradeModuleUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::actionButton(
    ns("action"), "Upgrade",
    width = "auto", class = "quick-button"
  )
}

UpgradeModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    showModal <- function() {
      body <- tagList(
        tags$iframe(
          src = "https://email.bigomics.ch/buy-now/", # Replace with the desired URL
          width = "100%",
          # height = "82vh",
          frameborder = "0"
        )
      )

      modal <- modalDialog2(
        title = NULL,
        # bsutils::modalHeader(
        #   div(class = "modal-title", "Share the Love. Invite A Friend."),
        #   style = "background-color: #f0f9fd;"
        # ),
        body,
        footer = NULL,
        size = "midscreen",
        easyClose = TRUE
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

    ## return
    list(
      click = ext_click ## exported function!
    )
  }) ## end of moduleServer
}
