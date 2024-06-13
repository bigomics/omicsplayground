##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

DatasetReportUI <- function(id) {
    ns <- shiny::NS(id)
    shiny::actionButton(
        ns("show_report_modal"),
        label = "Generate Report",
        icon = icon("file"),
        class = "btn btn-secondary",
        width = NULL
    )
}

DatasetReportServer <- function(
    id) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns ## NAMESPACE

        showModal <- function() {
            body <- tagList(
                div(
                    shiny::textInput(
                        inputId = ns("available_datasets"),
                        label = "", placeholder = "Select a dataset"
                    )
                ),
                shiny::actionButton(ns("generate_report_action"), "Submit", class = "btn btn-primary")
            )

            modal <- shiny::modalDialog(
                title = NULL,
                bsutils::modalHeader(
                    div(class = "modal-title", "Create a report")
                ),
                body,
                footer = NULL,
                size = "l",
                easyClose = TRUE
            )

            shiny::showModal(modal)
        }

        shiny::observeEvent(input$show_report_modal, {
            url <- Sys.getenv("API_BACKEND_URL")
            bearer_token <- Sys.getenv("BEARER_TOKEN")

            response <- httr::GET(
                glue::glue("{url}/secure/datasets"),
                httr::add_headers(Authorization = paste("Bearer", bearer_token))
            )

            print(httr::content(response, "text"))
            showModal()
        })
    }) ## end of moduleServer
}
