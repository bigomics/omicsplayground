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

    shiny::observeEvent(
      {
        input$action
      },
      {
        showModal()
      }
    )
  }) ## end of moduleServer
}
