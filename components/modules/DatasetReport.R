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
    auth,
    pgxtable) {
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

            dataset <- pgxtable$data()[pgxtable$rows_selected(), "dataset"]

            body <- tagList(
                div(
                    shiny::h4(paste("Generate report for dataset: ", dataset)),
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
                        choices = c("PDF", "HTML"),
                        selected = "PDF",
                        multiple = FALSE
                    ),
                    shiny::downloadButton(ns("download_pdf"), "Submit")
                )
            )

            output$download_pdf <- shiny::downloadHandler(
                filename = function() {
                    paste(dataset, paste0(".", input$output_format), sep = "")
                },
                content = function(file) {
                    shiny::removeModal()

                    shinyalert::shinyalert(
                        title = "Your report is being computed!",
                        text = "Your download will start shortly.",
                        type = "info"
                    )

                    # create a switch statement to replace pdf by poster-typst
                    output_format <- switch(input$output_format,
                        "PDF" = "poster-typst",
                        "HTML" = "html"
                    )

                    pgx_path <- file.path(auth$user_dir, paste0(dataset, ".pgx", ""))
                    print(pgx_path)

                    system2(
                        "quarto",
                        args = c(
                            "render",
                            file.path(quarto_file_path, "visreport.qmd"),
                            "--output",
                            "--to",
                            tolower(output_format),
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
                list(input$show_report_modal)
            },
            {
                ai <- 1

                req(pgxtable, !is.null(pgxtable$rows_selected()))
                dataset <- pgxtable$data()[pgxtable$rows_selected(), "dataset"]
                pgx <- playbase::pgx.load(file.path(auth$user_dir, paste0(dataset, ".pgx")))

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
