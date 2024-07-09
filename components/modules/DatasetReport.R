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
    id,
    auth) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns ## NAMESPACE

        quarto_file_path <- Sys.getenv("QUARTO_FILE_PATH")

        if (is.null(quarto_file_path)) {
            warning("QUARTO_FILE_PATH is not set. Please set it to the path where the quarto file is located.")
            return(NULL)
        }


        showModal <- function() {
            shiny::req(auth$logged)
            if (is.null(auth$logged) || !auth$logged) {
                return(NULL)
            }

            info <- playbase::pgxinfo.read(auth$user_dir, file = "datasets-info.csv")

            body <- tagList(
                div(
                    shiny::selectInput(
                        inputId = ns("available_datasets"),
                        label = "Select the dataset:",
                        choices = info$dataset,
                        selected = info$dataset[1]
                    ),
                    shiny::selectizeInput(
                        inputId = ns("sel_contrasts"),
                        label = "Select one or more comparisons",
                        choices = c("Getting comparisons..."),
                        selected = "Getting comparisons...",
                        multiple = TRUE
                    ),
                    shiny::selectizeInput(
                        inputId = ns("output_format"),
                        label = "Output format:",
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

                    pgx_path <- file.path(auth$user_dir, paste0(input$available_datasets, ".pgx", ""))
                    print(pgx_path)

                    system2(
                        "quarto",
                        args = c(
                            "render",
                            file.path(quarto_file_path, "visreport.qmd"),
                            "--output",
                            "--to",
                            tolower(input$output_format),
                            "-P",
                            paste("pgxdir:", pgx_path, sep = ""),
                            "-P",
                            paste("comparison:", input$sel_contrasts[1], sep = "")
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
        shiny::observeEvent(
            {
                list(input$available_datasets, input$show_report_modal)
            },
            {
                req(input$available_datasets)
                pgx <- playbase::pgx.load(file.path(auth$user_dir, paste0(input$available_datasets, ".pgx")))

                contrasts <- playbase::pgx.getContrasts(pgx)

                updateSelectizeInput(
                    session,
                    "sel_contrasts",
                    choices = contrasts,
                    selected = contrasts[1]
                )
            }
        ) # end of observe dataset and update contrasts

        # observe dataset and
    }) ## end of moduleServer
}
