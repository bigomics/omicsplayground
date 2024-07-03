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

            body <- tagList(
                div(
                    shiny::selectInput(
                        inputId = ns("available_datasets"),
                        label = "",
                        choices = datasets,
                        selected = datasets[1]
                    ),
                    shiny::selectizeInput(
                        inputId = ns("sel_contrasts"),
                        label = "Select one or more contrasts",
                        choices = c("Getting comparisons..."),
                        selected = "Getting comparisons...",
                        multiple = TRUE
                    ),
                    shiny::selectizeInput(
                        inputId = ns("output_format"),
                        label = "Output",
                        choices = c("PDF" = "pdf", "HTML" = "html"),
                        selected = "PDF",
                        multiple = FALSE
                    ),
                    shiny::downloadButton(ns("download_pdf"), "Submit")
                )
            )

            output$download_pdf <- shiny::downloadHandler(
                filename = function() {
                    paste(input$available_datasets, paste0(".", input$output_format), sep = "")
                },
                content = function(file) {
                    shiny::removeModal()

                    shinyalert::shinyalert(
                        title = "Your report is being computed!",
                        text = "Your download will start shortly.",
                        type = "info"
                    )

                    print("Generating report")
                    print(input$available_datasets)
                    print(input$sel_contrasts)

                    system2(
                        "quarto",
                        args = c(
                            "render",
                            file.path(quarto_file_path, "main.qmd"),
                            "--output",
                            "--to",
                            tolower(input$output_format),
                            "-P",
                            paste("dataset:", input$available_datasets, sep = ""),
                            "-P",
                            paste("comparisons:", paste0(input$sel_contrasts, collapse = ","), sep = "")
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

        ## observe dataset and update contrasts
        shiny::observeEvent(input$available_datasets, {
            print("Dataset changed")
            url <- Sys.getenv("API_BACKEND_URL")
            bearer_token <- Sys.getenv("BEARER_TOKEN")

            response <- httr::GET(
                glue::glue("{url}/secure/contrast_names?dataset={input$available_datasets}"),
                httr::add_headers(Authorization = paste("Bearer", bearer_token))
            )

            res <- jsonlite::fromJSON(httr::content(response, "text"))
            print(res)

            updateSelectizeInput(
                session,
                "sel_contrasts",
                choices = res,
                selected = res[1]
            )
        }) # end of observe dataset and update contrasts

        # observe dataset and
    }) ## end of moduleServer
}
