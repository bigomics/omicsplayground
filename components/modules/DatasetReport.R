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
            # Fetch the datasets from the API
            url <- Sys.getenv("API_BACKEND_URL")
            bearer_token <- Sys.getenv("BEARER_TOKEN")
            quarto_file_path <- Sys.getenv("QUARTO_FILE_PATH")

            response <- httr::GET(
                glue::glue("{url}/secure/datasets"),
                httr::add_headers(Authorization = paste("Bearer", bearer_token))
            )

            datasets <- jsonlite::fromJSON(httr::content(response, "text"))
            print(datasets)

            body <- tagList(
                div(
                    shiny::selectInput(
                        inputId = ns("available_datasets"),
                        label = "",
                        choices = datasets,
                        selected = datasets[1]
                    ),
                    shiny::downloadButton(ns("download_pdf"), "Submit")
                )
            )

            output$download_pdf <- shiny::downloadHandler(
                filename = function() {
                    paste(input$available_datasets, ".pdf", sep = "")
                },
                content = function(file) {
                    system2(
                        "quarto",
                        args = c(
                            "render",
                            file.path(quarto_file_path, "main.qmd"),
                            "--output",
                            "--to",
                            "pdf",
                            "-P",
                            paste("dataset:", input$available_datasets, sep = "")
                        ),
                        stdout = file
                    )
                }
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
            showModal()
        })
    }) ## end of moduleServer
}
