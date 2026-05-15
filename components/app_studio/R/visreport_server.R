##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

VisReportSettings <- function(id, output_format=NULL, type="dataset") {
  ns <- shiny::NS(id)

  info ="Create and download a summary report for the current dataset. The results in the report are created using default values and are intended only as for quick exploration of your data. For custom analysis, use the interactive dashboards."
  
  div.format <- NULL
  if(is.null(output_format)) {
    div.format <- shiny::tagList(
      br(),
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
      br()
    )
  }
  
  title <- ifelse(output_format=='poster', 'Poster', 'Slide deck')
    
  ui <- div(
    style = "padding: 10px 15px;",
    h3(title),    
    br(),
    div.format,
    shiny::selectizeInput(
      inputId = ns("report_type"),
      label = "Select report type:",
      choices = c(
        "Dataset summary" = "dataset",
        "Comparison summary" = "comparison"
      ),
      selected = type,
      multiple = FALSE
    ),
    shiny::br(),
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
    ),
    shiny::br(),
    shiny::h5("Select sections to include:"),
    shiny::checkboxInput(
      inputId = ns("section_clustering"),
      label = "Clustering (PCA, Heatmap)",
      value = TRUE
    ),
    shiny::conditionalPanel(
      "input.report_type == 'dataset'",
      ns = ns,
      shiny::checkboxInput(
        inputId = ns("section_phenotype"),
        label = "Phenotype Analysis",
        value = TRUE
      )
    ),
    shiny::checkboxInput(
      inputId = ns("section_differential"),
      label = "Differential Expression",
      value = TRUE
    ),
    shiny::checkboxInput(
      inputId = ns("section_enrichment"),
      label = "Enrichment Analysis",
      value = TRUE
    ),
    shiny::conditionalPanel(
      "input.report_type == 'comparison'",
      ns = ns,
      shiny::checkboxInput(
        inputId = ns("section_functional"),
        label = "Functional Analysis",
        value = TRUE
      ),
      shiny::checkboxInput(
        inputId = ns("section_drug"),
        label = "Drug Connectivity",
        value = TRUE
      )
    ),
    br(),
    br(),
    shiny::downloadButton(ns("download_pdf"), "Download",
      class = "btn btn-primary"
    )
  )
  return(ui)
}

VisReportUI <- function(id, output_format=NULL, type="dataset") {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("preview"))
}

VisReportServer <- function(id, pgx, output_format=NULL) {
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

    if (is.null(quarto_file_path) || quarto_file_path == "") {
      warning("[DatasetReportServer] ERROR: Please set QUARTO_FILE_PATH.")
      ## return(NULL)
    }

    r_format <- reactive({
      if(is.null(output_format)) {
        return(input$output_format)
      }
      return(output_format)
    })

    
    output$preview <- shiny::renderUI({
      format <- r_format()
      rtype <- input$report_type

      if (rtype == "dataset" && format == "poster") {
        thumb_img <- "visreport-dataset-poster-thumb.png"
      } else if (rtype == "dataset" && format == "slide") {
        thumb_img <- "visreport-dataset-slide-thumb.png"
      } else if (rtype == "comparison" && format == "poster") {
        thumb_img <- "visreport-comparison-poster-thumb.png"
      } else if (rtype == "comparison" && format == "slide") {
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
    
    output$download_pdf <- shiny::downloadHandler(
      filename = function() {
        ext <- switch(r_format(),
          "poster" = "pdf",
          "slide"  = "pdf",
          "report" = "html"
        )
        dataset <- pgx$name
        dataset <- sub("[.]pgx$", "", dataset)        
        paste0(dataset, "-", input$report_type, "-", r_format(), ".", ext)
      },
      content = function(file) {
        shiny::removeModal()

        if (is.null(quarto_file_path) || quarto_file_path == "") {
          shinyalert::shinyalert(
            title = "Duh!",
            text = "Reporting is not available for this server",
            type = "error",
            immediate = TRUE
          )
          return(NULL)
        }
        
        shinyalert::shinyalert(
          title = "Preparing your report!",
          text = "We are chopping the vegetables, boiling the water and simmering the sauce. Stay with us. Your download will start shortly.",
          type = "",
          imageUrl = "https://media4.giphy.com/media/v1.Y2lkPTc5MGI3NjExNTRpeG8wcTZyOXJxOTRramhsM3p3Z2wwNDVlbm5nNDdnZXlsc3ludyZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9cw/JSYv0MWRkzOWCk8f0Z/giphy.webp",
          imageWidth = 300,
          imageHeight = 200,
          closeOnEsc = FALSE
        )

        pgx_path <- dirname(pgx$filename)
        pgx_file <- basename(pgx$filename)
        sel_dataset <- sub("[.]pgx$","",pgx_file)        
        user <- pgx$creator

        ## create a switch statement to replace pdf by poster-typst
        render_format <- switch(r_format(),
          "poster"  = "poster-typst",
          "slide"   = "pdf",
          "report " = "html"
        )

        ## Create a Progress object
        progress <- shiny::Progress$new()

        ## quarto does not like output with rel/abs
        ## path. Also switching to temp folder is safer to
        ## avoid multiple users using quarto at the same time.
        cur_wd <- getwd()
        
        on.exit({
          progress$close()
          setwd(cur_wd)
        })
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
          sel_contrasts <- playbase::pgx.getContrasts(pgx)
        }
        
        progress$set(message = "Creating report", value = 0)

        if (input$report_type %in% c("dataset", "report")) {
          progress$inc(0.3, detail = "Creating dataset summary...")

          qmd_file <- "visreport-dataset.qmd"
          if (r_format() == "slide") {
            qmd_file <- "visreport-bigpage.qmd"
          }
          all_ct <- paste(sel_contrasts, collapse = ",")
          
          # Prepare section parameters for dataset reports
          section_args <- c(
            "-P", paste0("sections_clustering:", tolower(input$section_clustering)),
            "-P", paste0("sections_phenotype:", tolower(input$section_phenotype)),
            "-P", paste0("sections_differential:", tolower(input$section_differential)),
            "-P", paste0("sections_enrichment:", tolower(input$section_enrichment))
          )
          
          system2(
            "quarto",
            args = c(
              "render",
              file.path(quarto_file_path, qmd_file),
              "--output -",
              paste("--to poster-typst"),
              "-P", shQuote(paste0("pgxdir:", pgx_path)),
              "-P", shQuote(paste0("comparisons:", all_ct)),
              "-P", shQuote(paste0("dataset:", sel_dataset)),
              "-M", shQuote(paste0("user:", user)),
              "-M", shQuote(paste0("title:", sel_dataset)),
              section_args
            ),
            stdout = file
          )
          
        }

        if (input$report_type == "comparison") {
          files <- c()
          ncontrasts <- length(sel_contrasts)
          for (i in 1:ncontrasts) {
            ct <- sel_contrasts[i]
            progress$inc(
              1 / (ncontrasts + 1),
              detail = paste0("summarizing comparison ", i, "/", ncontrasts)
            )
            tmp <- tempfile()

            # Prepare section parameters for comparison report
            section_args <- c(
              "-P", paste0("sections_clustering:", tolower(input$section_clustering)),
              "-P", paste0("sections_differential:", tolower(input$section_differential)),
              "-P", paste0("sections_enrichment:", tolower(input$section_enrichment)),
              "-P", paste0("sections_functional:", tolower(input$section_functional)),
              "-P", paste0("sections_drug:", tolower(input$section_drug))
            )

            system2(
              "quarto",
              args = c(
                "render",
                file.path(quarto_file_path, "visreport-comparison.qmd"),
                "--output -",
                paste("--to", render_format),
                "-P", shQuote(paste0("pgxdir:", pgx_path)),
                "-P", shQuote(paste0("comparison:", ct)),
                "-P", shQuote(paste0("dataset:", sel_dataset)),
                "-M", shQuote(paste0("dataset:", sel_dataset)),
                "-M", shQuote(paste0("user:", user)),
                "-M", shQuote(paste0("title:", ct)),
                section_args
              ),
              stdout = tmp
            )
            files <- c(files, tmp)
          }
          ## remove any empty files
          zero_size_files <- files[file.size(files) == 0]
          if (length(zero_size_files) > 0) {
            files <- files[file.size(files) > 0]
            message("[DatasetReportServer:download_pdf] removed ", length(zero_size_files), " empty file(s)")
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

        record_report_download(input$report_type)
        
        shinyalert::shinyalert(
          title = "Ready!",
          text = paste("Your",r_format(),"is finished. Please check your downloads folder."),
          type = "",
          immediate = TRUE
        )

        ## for(i in 1:length(files)) unlink(files[i])
      } ## end-of-content
    )  ## end of download_pdf

    ## observe dataset and update contrasts
    shiny::observeEvent({
      list(pgx$X)
    }, {
      contrasts <- playbase::pgx.getContrasts(pgx)
      shiny::req(contrasts)
      updateSelectizeInput(
        session,
        "sel_contrasts",
        choices = contrasts,
        selected = head(contrasts, 6)
      )
    }) # end of observe dataset and update contrasts

    # observe dataset and
  }) ## end of moduleServer
}
