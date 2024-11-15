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
    if (quarto_file_path == "") {
      ## default to parent of OPG
      path1 <- file.path(OPG, "../pgx-visreport")
      path2 <- file.path(OPG, "./libx/pgx-visreport")
      if (dir.exists(path1)) quarto_file_path <- path1
      if (dir.exists(path2)) quarto_file_path <- path2
    }

    dbg("[DatasetReportServer] quarto_file_path = ", quarto_file_path)

    if (is.null(quarto_file_path) || quarto_file_path == "") {
      warning("[DatasetReportServer] ERROR: Please set QUARTO_FILE_PATH.")
      ## return(NULL)
    }

    output$visreport_thumbnail <- shiny::renderUI({
      format <- input$output_format
      rtype <- input$report_type

      if (rtype == "dataset-summary" && format == "poster") {
        thumb_img <- "visreport-dataset-poster-thumb.png"
      } else if (rtype == "dataset-summary" && format == "slide") {
        thumb_img <- "visreport-dataset-slide-thumb.png"
      } else if (rtype == "comparison-summary" && format == "poster") {
        thumb_img <- "visreport-comparison-poster-thumb.png"
      } else if (rtype == "comparison-summary" && format == "slide") {
        thumb_img <- "visreport-comparison-slide-thumb.png"
      } else {
        ## never here...
        thumb_img <- "monster-hi.png"
      }

      shiny::tags$img(
        src = paste0("static/", thumb_img),
        width = "100%"
      )
    })

    showModal <- function() {
      shiny::req(auth$logged)
      if (is.null(auth$logged) || !auth$logged) {
        return(NULL)
      }

      dataset <- pgxtable$data()[pgxtable$rows_selected(), "dataset"]
      dataset <- sub("[.]pgx$", "", dataset)

      if (length(dataset) == 0 || is.null(dataset) || dataset == "") {
        shinyalert::shinyalert(
          title = "",
          text = "Please select a dataset from the table first",
          type = "error",
          closeOnEsc = TRUE
        )
        return(NULL)
      }

      body <- shiny::tagList(
        "Create and download a summary report for the current dataset '",
        shiny::tags$b(dataset), "'.",
        "The results in the report are created using default values and are intended only as for quick exploration of your data. For custom analysis, use the interactive dashboards.",
        shiny::br(), shiny::br(),
        shiny::hr(), shiny::br(),
        bslib::layout_columns(
          col_widths = c(4, 8),
          shiny::div(
            shiny::selectizeInput(
              inputId = ns("report_type"),
              label = "Select report type:",
              choices = c(
                "Dataset summary" = "dataset-summary",
                "Comparison summary" = "comparison-summary"
              ),
              selected = 1,
              multiple = FALSE
            ),
            shiny::br(),
            shiny::selectizeInput(
              inputId = ns("output_format"),
              label = "Select output format:",
              choices = c(
                "Poster (PDF, one big page)" = "poster",
                "Slides (PDF, multiple pages)" = "slide"
                ## "Document (HTML)" = "report"
              ),
              selected = 1,
              multiple = FALSE
            ),
            br(),
            shiny::radioButtons(
              inputId = ns("contrasts_choice"),
              label = "Select comparisons to include in your report:",
              choices = c("all (maximum 20)" = "all", "select")
            ),
            shiny::conditionalPanel(
              "input.contrasts_choice == 'select'",
              ns = ns,
              shiny::selectizeInput(
                inputId = ns("sel_contrasts"),
                label = "Select comparisons:",
                choices = c("Getting comparisons..."),
                selected = "Getting comparisons...",
                multiple = TRUE
              )
            )
          ),
          div(
            shiny::uiOutput(ns("visreport_thumbnail"))
          )
        )
      )

      output$download_pdf <- shiny::downloadHandler(
        filename = function() {
          ext <- switch(input$output_format,
            "poster" = "pdf",
            "slide"  = "pdf",
            "report" = "html"
          )
          paste0(dataset, "-", input$report_type, "-", input$output_format, ".", ext)
        },
        content = function(file) {
          shiny::removeModal()

          shinyalert::shinyalert(
            title = "Preparing your report!",
            text = "We are chopping the vegetables, boiling the water and simmering the sauce. Stay with us. Your download will start shortly.",
            type = "",
            imageUrl = "https://media4.giphy.com/media/v1.Y2lkPTc5MGI3NjExNTRpeG8wcTZyOXJxOTRramhsM3p3Z2wwNDVlbm5nNDdnZXlsc3ludyZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9cw/JSYv0MWRkzOWCk8f0Z/giphy.webp",
            imageWidth = 300,
            imageHeight = 200,
            closeOnEsc = FALSE
          )

          pgx_path <- auth$user_dir
          sel_dataset <- pgxtable$data()[pgxtable$rows_selected(), "dataset"]
          sel_dataset <- sub("[.]pgx$", "", sel_dataset)
          pgx_file <- paste0(sel_dataset, ".pgx")

          user <- ifelse(shiny::isTruthy(auth$email), auth$email, "-")

          ## create a switch statement to replace pdf by poster-typst
          render_format <- switch(input$output_format,
            "poster"  = "poster-typst",
            "slide"   = "pdf",
            "report " = "html"
          )

          ## Create a Progress object
          progress <- shiny::Progress$new()
          on.exit(progress$close())

          ## quarto does not like output with rel/abs
          ## path. Also switching to temp folder is safer to
          ## avoid multiple users using quarto at the same time.
          cur_wd <- getwd()
          qdir <- quarto_file_path
          tmpdir <- tempdir()
          setwd(tmpdir)
          file.copy(file.path(qdir, "visreport-dataset.qmd"), ".")
          file.copy(file.path(qdir, "visreport-comparison.qmd"), ".")
          file.copy(file.path(qdir, "visreport-bigpage.qmd"), ".")
          file.copy(file.path(qdir, "_extensions"), ".", recursive = TRUE)
          file.copy(file.path(qdir, "images"), ".", recursive = TRUE)

          sel_contrasts <- input$sel_contrasts
          if (input$contrasts_choice == "all") {
            sel_contrasts <- getContrasts()
          }
          dbg("[DatasetReportServer] sel_contrasts = ", sel_contrasts)

          progress$set(message = "Creating report", value = 0)

          if (input$report_type %in% c("dataset-summary", "report")) {
            progress$inc(0.3, detail = "Creating dataset summary...")

            qmd_file <- "visreport-dataset.qmd"
            if (input$output_format == "slide") {
              qmd_file <- "visreport-bigpage.qmd"
            }
            all_ct <- paste(sel_contrasts, collapse = ",")
            system2(
              "quarto",
              args = c(
                "render",
                file.path(quarto_file_path, qmd_file),
                "--output -",
                paste("--to poster-typst"),
                "-P", paste0("pgxdir:", pgx_path),
                "-P", paste0("comparisons:", all_ct),
                "-P", paste0("dataset:", sel_dataset),
                "-M", paste0("user:", user),
                "-M", paste0("title:", sel_dataset)
              ),
              stdout = file
            )
          }

          if (input$report_type == "comparison-summary") {
            files <- c()
            ncontrasts <- length(sel_contrasts)
            for (i in 1:ncontrasts) {
              ct <- sel_contrasts[i]
              progress$inc(
                1 / (ncontrasts + 1),
                detail = paste0("summarizing comparison ", i, "/", ncontrasts)
              )
              tmp <- tempfile()
              system2(
                "quarto",
                args = c(
                  "render",
                  file.path(quarto_file_path, "visreport-comparison.qmd"),
                  "--output -",
                  paste("--to", render_format),
                  "-P", paste0("pgxdir:", pgx_path),
                  "-P", paste0("comparison:", ct),
                  "-P", paste0("dataset:", sel_dataset),
                  "-M", paste0("dataset:", sel_dataset),
                  "-M", paste0("user:", user),
                  "-M", paste0("title:", ct)
                ),
                stdout = tmp
              )
              files <- c(files, tmp)
            }

            ## finally merge all pages
            progress$inc(1 / (ncontrasts + 1), detail = "Merging pages")
            message("[DatasetReportServer:download_pdf] merging PDF pages...")
            system2(
              "pdftk",
              paste(paste(files, collapse = " "), "cat output", file)
            )
            for (i in 1:length(files)) unlink(files[i])
          }

          ## return to app running folder
          setwd(cur_wd)
          unlink(tmpdir)

          shinyalert::shinyalert(
            title = "Yay! Your report is ready",
            text = "We finished your report. Please check your downloads folder.",
            type = "",
            immediate = TRUE
          )

          ## for(i in 1:length(files)) unlink(files[i])
        } ## end-of-content
      )

      modal <- shiny::modalDialog(
        title = NULL,
        bsutils::modalHeader(
          div(class = "modal-title", "Create Summary Report")
        ),
        body,
        ## footer = NULL,
        footer = shiny::tagList(
          shiny::modalButton("Cancel"),
          shiny::downloadButton(ns("download_pdf"), "Download",
            class = "btn btn-primary"
          )
        ),
        size = "xl",
        easyClose = FALSE
      )

      shiny::showModal(modal)
    }

    shiny::observeEvent(input$show_report_modal, {
      if (is.null(quarto_file_path) || quarto_file_path == "") {
        shinyalert::shinyalert(
          title = "Duh!",
          text = "Reporting is not available for this server",
          type = "error",
          immediate = TRUE
        )
      } else {
        showModal()
      }
    })

    getContrasts <- eventReactive(
      {
        list(input$show_report_modal)
      },
      {
        req(pgxtable, !is.null(pgxtable$rows_selected()))
        dataset <- pgxtable$data()[pgxtable$rows_selected(), "dataset"]
        pgx <- playbase::pgx.load(file.path(auth$user_dir, paste0(dataset, ".pgx")))
        playbase::pgx.getContrasts(pgx)
      }
    )

    ## observe dataset and update contrasts
    shiny::observeEvent(
      {
        list(input$show_report_modal)
      },
      {
        contrasts <- getContrasts()
        shiny::req(contrasts)

        updateSelectizeInput(
          session,
          "sel_contrasts",
          choices = contrasts,
          selected = head(contrasts, 6)
        )
      }
    ) # end of observe dataset and update contrasts

    # observe dataset and
  }) ## end of moduleServer
}
