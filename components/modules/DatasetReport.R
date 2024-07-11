##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

DatasetReportUI <- function(id) {
    ns <- shiny::NS(id)
    shiny::actionButton(
        ns("show_report_modal"),
        label = "Generate report",
        icon = icon("file"),
        class = "btn btn-outline-secondary",
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
                    "Download a summary report for the current dataset ", shiny::tags$b(dataset),
                    shiny::br(),shiny::br(),
                    shiny::selectizeInput(
                        inputId = ns("report_type"),
                        label = "Select report type:",
                        choices = c(
                          "Dataset summary" = "dataset-summary",
                          "Comparison summary" = "comparison-summary",
                        ),
                        selected = 1,
                        multiple = FALSE
                    ),
                    shiny::selectizeInput(
                        inputId = ns("output_type"),
                        label = "Select output format:",
                        choices = c(
                          "Poster (PDF, one big page)" = "poster",
                          "Slides (PDF, multiple pages)" = "slides",
                          "Document (HTML)" = "html"
                        ),
                        selected = 1,
                        multiple = FALSE
                    ),
                    br(),
                    shiny::selectizeInput(
                        inputId = ns("sel_contrasts"),
                        label = "Select comparisons to include in your report:",
                        choices = c("Getting comparisons..."),
                        selected = "Getting comparisons...",
                        multiple = TRUE
                      )
                    ##shiny::downloadButton(ns("download_pdf"), "Submit")
                )
            )

            output$download_pdf <- shiny::downloadHandler(
                filename = function() {
                    ext <- switch(
                        input$report_type,
                        "dataset-summary"  = "pdf",
                        "comparison-summary"  = "pdf",
                        "html-report" = "html"
                    )                                              
                    paste0(dataset, "-", input$report_type, ".", ext)
                },
                content = function(file) {
                    shiny::removeModal()

                    shinyalert::shinyalert(
                        title = "Preparing your report!",
                        text = "We are chopping the vegetables, boiling the water and simmering the sauce. Stay with us. Your download will start shortly.",
                        type = "info",
                        closeOnEsc = FALSE
                    )
                    
                    pgx_path <- auth$user_dir
                    sel_dataset <- pgxtable$data()[pgxtable$rows_selected(), "dataset"]
                    sel_dataset <- sub("[.]pgx$","",sel_dataset)
                    pgx_file <- paste0(sel_dataset, ".pgx")

                    user <- ifelse( shiny::isTruthy(auth$email), auth$email, "-")
                    dbg("[DatasetReportServer:download_pdf] user = ",user)
                    
                    ## create a switch statement to replace pdf by poster-typst
                    render_format <- switch(
                        input$report_type,
                        "dataset-summary"  = "poster-typst",
                        "comparison-summary"  = "poster-typst",
                        "html-report" = "html"
                    )

                    ## Create a Progress object
                    progress <- shiny::Progress$new()
                    on.exit(progress$close())
                    
                    ## quarto does not like output with rel/abs
                    ## path. Also switching to temp folder is safer to
                    ## avoid multiple users using quarto at the same time.
                    cur_wd <- getwd()
                    qdir <- quarto_file_path
                    tmpdir = tempdir()
                    setwd(tmpdir)
                    file.copy( file.path(qdir,"visreport-dataset.qmd"), "." )
                    file.copy( file.path(qdir,"visreport-comparison.qmd"), "." )
                    file.copy( file.path(qdir,"_extensions"), ".", recursive = TRUE )
                    file.copy( file.path(qdir,"images"), ".", recursive = TRUE )
                    
                    progress$set(message = "Creating report", value = 0)

                    if(input$report_type == "dataset-summary") {
                      progress$inc(0.3, detail = "Creating dataset summary...")
                      all_ct <- paste(input$sel_contrasts, collapse=",")
                      system2(
                        "quarto",
                        args = c(
                          "render",
                          file.path(quarto_file_path, "visreport-dataset.qmd"),
                          "--output -",
                          "--to poster-typst",
                          "-P", paste0("pgxdir:", pgx_path),
                          "-P", paste0("comparisons:", all_ct),
                          "-P", paste0("dataset:", sel_dataset),
                          "-M", paste0("user:", user),
                            "-M", paste0("title:",sel_dataset)
                        ),
                        stdout = file
                      )
                    }

                    if(input$report_type == "comparison-summary") {                    
                      files <- c()                    
                      ncontrasts <- length(input$sel_contrasts)
                      for(i in 1:ncontrasts) {
                        ct <- input$sel_contrasts[i]
                        progress$inc(
                           1/(ncontrasts+1),
                           detail = paste0("summarizing comparison ",i,"/",ncontrasts)
                        )
                        tmp <- tempfile()
                        system2(
                          "quarto",
                          args = c(
                            "render",
                            file.path(quarto_file_path, "visreport-comparison.qmd"),
                            "--output -",
                            "--to poster-typst",
                            "-P", paste0("pgxdir:", pgx_path),
                            "-P", paste0("comparison:", ct),
                            "-P", paste0("dataset:", sel_dataset),
                            "-M", paste0("user:", user),
                            "-M", paste0("title:", ct)                          
                          ),
                        stdout = tmp
                        )
                        files <- c(files, tmp)                      
                      }

                      ## finally merge all pages
                      progress$inc(1/(ncontrasts+1), detail = "Merging pages")
                      message("[DatasetReportServer:download_pdf] merging PDF pages...")
                      system2(
                        "pdftk",
                        paste(paste(files,collapse=" "),"cat output",file)
                      )
                      for(i in 1:length(files)) unlink(files[i])
                    }

                    if(input$report_type == "html-report") {
                      progress$inc(0.3, detail = "Creating HTML report...")
                      all_ct <- paste(input$sel_contrasts, collapse=",")
                      system2(
                        "quarto",
                        args = c(
                          "render",
                          "visreport-dataset.qmd",
                          paste("--output report.html"),
                          "--to html",
                          "-P", paste0("pgxdir:", pgx_path),
                          "-P", paste0("comparisons:", all_ct),
                          "-P", paste0("dataset:", sel_dataset),
##                          "-P", paste0("user:", user),
                          "-M", paste0("user:", user),
                          "-M", paste0("title:",sel_dataset)
                        )
                      )
                      file.copy("report.html", file)
                    }

                    ## return to app running folder
                    setwd(cur_wd)
                    unlink(tmpdir)

                    dbg("[DatasetReportServer:download_pdf] done!")
                    shinyalert::shinyalert(
                        title = "Yay! Your report is ready",
                        text = "We finished creating your report. Please check your downloads folder.",
                        type = "info",
                        immediate = TRUE
                    )                    

                    ## for(i in 1:length(files)) unlink(files[i])
                } ## end-of-content
                
            )

            modal <- shiny::modalDialog(
                title = NULL,
                bsutils::modalHeader(
                    div(class = "modal-title", "Create summary report")
                ),
                body,
                ## footer = NULL,
                footer = shiny::tagList(
                   shiny::modalButton("Cancel"),
                   shiny::downloadButton(ns("download_pdf"), "Download", class=NULL)
                ),
                size = "l",
                easyClose = FALSE
            )

            shiny::showModal(modal)
        }

        shiny::observeEvent(input$show_report_modal, {
            showModal()
        })

        ## observe dataset and update contrasts
        shiny::observeEvent({
            list(input$show_report_modal)
        }, {
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
        }) # end of observe dataset and update contrasts

        # observe dataset and
    }) ## end of moduleServer
}
