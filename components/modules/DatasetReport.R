##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

DatasetReportUI <- function(id) {
    ns <- shiny::NS(id)
    shiny::actionButton(
        ns("show_report_modal"), "Generate Report",
        width = "auto", class = "quick-button"
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
                # body,
                footer = NULL,
                size = "l",
                easyClose = TRUE,
                tags$style(".modal-dialog {width: 720px;}"),
                tags$style(".modal-content {background-color: #f0f9fd;}"),
                tags$style(".modal-header {padding: 0px;}")
            )

            shiny::showModal(modal)
        }

        shiny::observeEvent(input$show_report_modal, {
            print("generate report clicked")
            showModal()
        })
    }) ## end of moduleServer
}
