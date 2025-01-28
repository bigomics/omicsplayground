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

UpgradeModuleServer <- function(id, auth) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    showModal <- function() {
      body <- tagList(
        tags$iframe(
          src = "https://upgrade.bigomics.ch/buy-now/", # Replace with the desired URL
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

    allowed_domains <- shiny::reactive({
      domains_path <- paste0(ETC, "/allowed_domains.csv")
      if (file.exists(domains_path)) {
        logged_in_domain <- strsplit(auth$email, "@")[[1]][2]
        if (is.na(logged_in_domain)) logged_in_domain <- ""
        domain_found <- length(system(paste0("grep '", logged_in_domain, "' ", domains_path), intern = TRUE)) > 0
        return(domain_found)
      } else {
        return(FALSE)
      }
    })

    shiny::observeEvent(
      {
        input$action
      },
      {
        record_UPGRADE()
        if (allowed_domains()) {
          showModal()
        } else {
          browseURL("https://bigomics.ch/pricing/")
        }
      }
    )
  }) ## end of moduleServer
}
