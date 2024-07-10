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
        if(quarto_file_path == "") {
          quarto_file_path <- file.path(OPG, "../pgx-visreport")
        }

        if (is.null(quarto_file_path) || quarto_file_path == "") {
            warning("[DatasetReportServer] ERROR: Please set QUARTO_FILE_PATH to quarto file location.")
            return(NULL)
        }

        showModal <- function() {
            shiny::req(auth$logged)
            if (is.null(auth$logged) || !auth$logged) {
                return(NULL)
            }

            dataset <- pgxtable$data()[pgxtable$rows_selected(), "dataset"]
            dataset <- sub("[.]pgx$","",dataset)

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
                        choices = c(
                          "Visual summary (PDF)" = "pdf",
                          "HTML report" = "html"
                        ),
                        selected = "PDF",
                        multiple = FALSE
                    ),
                    shiny::downloadButton(ns("download_pdf"), "Submit")
                )
            )

            output$download_pdf <- shiny::downloadHandler(
                filename = function() {
                    paste0(dataset, "-report.", input$output_format)
                },
                content = function(file) {
                    shiny::removeModal()

                    shinyalert::shinyalert(
                        title = "Your report is being created!",
                        text = "Please wait. Your download will start shortly.",
                        type = "info"
                    )

                    
                    pgx_path <- auth$user_dir
                    sel_dataset <- pgxtable$data()[pgxtable$rows_selected(), "dataset"]
                    sel_dataset <- sub("[.]pgx$","",sel_dataset)
                    pgx_file <- paste0(sel_dataset, ".pgx")
                    
                    dbg("pgx_path = ",pgx_path)
                    dbg("pgx_file = ",pgx_file)
                    dbg("sel_dataset = ",sel_dataset)
                    dbg("names(auth) = ",names(auth))
                    dbg("auth$email = ",auth$email)
                    dbg("auth$username = ",auth$username)
                    
                    ## create a switch statement to replace pdf by poster-typst
                    render_format <- switch(input$output_format,
                        "pdf" = "poster-typst",
                        "html" = "html"
                    )

                    ## Create a Progress object
                    progress <- shiny::Progress$new()
                    on.exit(progress$close())
                    
                    progress$set(message = "Creating dataset summary...", value = 0)
                    all_ct <- paste(input$sel_contrasts, collapse=",")                    
                    files <- c()
                    tmp <- tempfile(fileext=".pdf")  ## ??                   
                    system2(
                      "quarto",
                      args = c(
                        "render",
                        file.path(quarto_file_path, "visreport-dataset.qmd"),
                        "--output -",
                        "--to", render_format,
                        "-P", paste0("pgxdir:", pgx_path),
                        "-P", paste0("comparison:", all_ct),
                        "-P", paste0("dataset:", pgx_file),
                        "-P", paste0("user:", auth$email)
                      ),
                      stdout = tmp
                    )
                    files <- c(files, tmp)
                    
                    ncontrasts <- length(input$sel_contrasts)
                    for(i in 1:ncontrasts) {
                      ct <- input$sel_contrasts[i]
                      progress$inc(1/(ncontrasts+1),
                                   detail = paste0("Creating summary for comparison ",i,"/",ncontrasts))

                      dbg("pgx_path = ",pgx_path)
                      dbg("pgx_file = ",pgx_file)
                      dbg("ct = ",ct)
                      dbg("auth$email = ",auth$email)

                      tmp <- tempfile(fileext=".pdf")
                      system2(
                        "quarto",
                        args = c(
                          "render",
                          file.path(quarto_file_path, "visreport-comparison.qmd"),
                          "--output -",
                          "--to", render_format,
                          "-P", paste0("pgxdir:", pgx_path),
                          "-P", paste0("comparison:", ct),
                          "-P", paste0("dataset:", pgx_file),
                          "-P", paste0("user:", auth$email)
                        ),
                        stdout = tmp
                      )
                      files <- c(files, tmp)                      
                    }

                    ## finally merge all pages
                    progress$inc(1/(ncontrasts+1), detail = "Merging pages")                      
                    if(render_format == "poster-typst") {
                      dbg("[DatasetReportServer:download_pdf] merging PDF pages...")
                      system2(
                        "pdftk",
                        paste(paste(files,collapse=" "),"cat output",file)
                      )
                    }
                    if(render_format == "html") {
                      dbg("[DatasetReportServer:download_pdf] merging HTML pages...")
                      system2(
                        "cat",
                        paste(files, collapse=" ")
                      )
                      stdout = file
                    }                    

                    dbg("[DatasetReportServer:download_pdf] done!")
                    shinyalert::shinyalert(
                        title = "Your report is ready!",
                        text = "Your report has been downloaded.",
                        type = "info",
                        immediate = TRUE
                    )                    
                    
                } ## end-of-content


                
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
                    selected = head(contrasts,6)
                )
            }
        ) # end of observe dataset and update contrasts

        # observe dataset and
    }) ## end of moduleServer
}
