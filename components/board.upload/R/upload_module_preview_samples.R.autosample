##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

upload_table_preview_samples_ui <- function(id) {
  ns <- shiny::NS(id)
  uiOutput(ns("table_samples"), fill = TRUE)
}

upload_table_preview_samples_server <- function(
    id,
    uploaded,
    checklist,
    scrollY,
    width,
    height,
    title,
    info.text,
    caption) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    table_data <- shiny::reactive({
      shiny::req(!is.null(uploaded$samples.csv))
      dt <- uploaded$samples.csv
      nrow0 <- nrow(dt)
      ncol0 <- ncol(dt)
      MAXSHOW <- 100
      if (nrow(dt) > MAXSHOW) {
        dt <- head(dt, MAXSHOW)
        dt <- rbind(dt, rep(NA, ncol(dt)))
        n1 <- nrow0 - MAXSHOW
        rownames(dt)[nrow(dt)] <- paste0("[+", n1, " rows]")
      }
      if (ncol(dt) > MAXSHOW) {
        dt <- dt[, 1:MAXSHOW]
        dt <- cbind(dt, rep(NA, nrow(dt)))
        n1 <- ncol0 - MAXSHOW
        colnames(dt)[ncol(dt)] <- paste0("[+", n1, " columns]")
      }
      dt
    })

    table.RENDER <- function() {
      dt <- table_data()
      req(!is.null(dt))

      DT::datatable(dt,
        class = "compact hover",
        rownames = TRUE,
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "column", selected = 1),
        options = list(
          dom = "lfrtip",
          scroller = TRUE,
          scrollX = TRUE,
          scrollY = scrollY,
          scrollResize = TRUE,
          deferRender = TRUE
        )
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    output$table_samples <- shiny::renderUI({
      action_buttons <- div(
        style = "display: flex; justify-content: left; margin-bottom: 8px;",
        shiny::actionButton(
          ns("create_from_samplenames"),
          "Create from sample names",
          class = "btn-sm btn-outline-primary m-1"
        ),
        div(
          if (!is.null(uploaded$samples.csv)) {
            shiny::actionButton(
              ns("remove_samples"),
              "Cancel",
              icon = icon("trash-can"),
              class = "btn-sm btn-outline-danger m-1"
            )
          } else {
            shiny::actionButton(
              ns("load_example"), "Load example data",
              class = "btn-sm btn-outline-primary m-1"
            )
          }
        ),
        shiny::actionButton(
          ns("check_documentation_samples"),
          "Read documentation",
          class = "btn-sm btn-outline-primary m-1",
          onclick = "window.open('https://omicsplayground.readthedocs.io/en/latest/dataprep/samples/', '_blank')"
        )
      )

      div(
        bslib::as_fill_carrier(),
        if (is.null(uploaded$samples.csv)) {
          bslib::layout_columns(
            col_widths = c(-3, 6, -3),
            row_heights = list("auto", 11, 1),
            gap = "0.3rem",
            bslib::as_fill_carrier(
              bs_alert("The samples file (samples.csv) contains the phenotypic information of your samples. The file should be a tabular text file (csv) with the samples in the rows and the phenotypic data (metadata) in the columns. The first column contains the sample names, which must be unique, and has to match the names given in the header of the counts file.", closable = FALSE),
              style = "align-items: end"
            ),
            bslib::card(
              fileInputArea(
                ns("samples_csv"),
                shiny::h4("Choose samples.csv", class = "mb-0"),
                multiple = FALSE,
                accept = c(".csv"),
                width = "100%"
              ),
              style = "background-color: aliceblue; border: 0.07rem dashed steelblue;"
            ),
            action_buttons
          )
        } else {
          bslib::layout_columns(
            col_widths = 12,
            bslib::layout_columns(
              col_widths = c(8, 4),
              TableModuleUI(
                ns("samples_datasets"),
                width = width,
                height = height,
                title = title,
                info.text = info.text,
                caption = caption,
                label = "",
                show.maximize = FALSE
              ),
              bslib::card(
                bslib::navset_pill(
                  bslib::nav_panel(
                    title = "UMAP",
                    br(),
                    plotOutput(ns("umap"), height = "500px")
                  ),
                  bslib::nav_panel(
                    title = "Distribution",
                    br(),
                    plotOutput(ns("phenotype_stats"), height = "500px")
                  )
                )
              )
            ),
            bslib::layout_columns(
              action_buttons,
              br(),
              uiOutput(ns("error_summary"))
            )
          )
        }
      )
    })

    output$error_summary <- renderUI({
      chk1 <- checklist$samples.csv$checks
      chk2 <- checklist$samples_counts$checks

      div(
        style = "display: flex; justify-content: right; vertical-align: text-bottom; margin: 8px;",
        check_to_html(
          checklist$samples.csv$checks,
          ## pass_msg = "All samples checks passed",
          pass_msg = "Samples matrix OK. ",
          null_msg = "Samples checks not run yet.
                            Fix any errors with samples first. ",
          false_msg = "Samples checks: warning ",
          details = FALSE
        ),
        check_to_html(
          checklist$samples_counts$checks,
          ##          pass_msg = tspan("All samples-counts checks passed"),
          pass_msg = tspan("Samples-counts OK.", js = FALSE),
          null_msg = tspan("Samples-counts checks not run yet.
                        Fix any errors with samples or counts first.", js = FALSE),
          false_msg = tspan("Samples-counts checks: warning", js = FALSE),
          details = FALSE
        )
      )
    })

    output$phenotype_stats <- renderPlot({
      pheno <- uploaded$samples.csv
      shiny::req(nrow(pheno))
      plotPhenoDistribution(data.frame(pheno))
    })

    output$umap <- renderPlot({
      counts <- uploaded$counts.csv
      shiny::req(nrow(counts))
      counts <- playbase::pgx.countNormalization(counts, "median.center.nz")
      prior <- min(counts[which(counts > 0)], na.rm = TRUE)
      X <- log2(counts + prior)
      Y <- uploaded$samples.csv
      sel <- grep("group|condition", colnames(Y), ignore.case = TRUE)
      sel <- head(c(sel, 1), 1)
      y <- Y[, sel]
      hilight2 <- colnames(X)
      if (ncol(X) > 100) hilight2 <- NULL
      shiny::validate(shiny::need(
        any(colnames(X) %in% rownames(Y)),
        "No matches between samples and counts."
      ))
      playbase::pgx.dimPlot(
        X, y,
        method = "umap",
        plotlib = "base",
        cex = 2.5,
        xlab = "umap-x",
        ylab = "umap-y",
        hilight2 = hilight2 ## label all points
      )
    })

    # error pop-up alert
    observeEvent(
      {
        list(
          checklist$samples.csv$checks,
          checklist$samples_counts$checks
        )
      },
      {
        checks1 <- checklist$samples.csv$checks
        checks2 <- checklist$samples_counts$checks

        if (length(checks1) > 0 || length(checks2) > 0) {
          err.html <- ""
          if (length(checks1) > 0) {
            err1 <- check_to_html(
              checks1,
              pass_msg = tspan("All samples checks passed", js = FALSE),
              null_msg = tspan("Samples checks not run yet.
                            Fix any errors with counts first.", js = FALSE),
              details = TRUE
            )
            err.html <- paste(err.html, err1)
          }
          if (length(checks2) > 0) {
            err2 <- check_to_html(
              checks2,
              pass_msg = tspan("All samples-counts checks passed", js = FALSE),
              null_msg = tspan("Samples-counts checks not run yet.
                        Fix any errors with samples or counts first.", js = FALSE),
              details = TRUE
            )
            err.html <- paste(err.html, err2)
          }
          shinyalert::shinyalert(
            title = "Warning",
            text = err.html,
            html = TRUE
          )
        }
      }
    )

    ## pass counts to uploaded when uploaded
    observeEvent(input$samples_csv, {
      # check if samples is csv (necessary due to drag and drop of any file)
      ext <- tools::file_ext(input$samples_csv$name)[1]
      if (ext != "csv") {
        shinyalert::shinyalert(
          title = "File format not supported.",
          text = "Please make sure the file is a CSV file.",
          type = "error"
        )
        return()
      }
      # if samples not in file name, give warning and return
      if (!grepl("sample", input$samples_csv$name, ignore.case = TRUE)) {
        shinyalert::shinyalert(
          title = "Samples not in filename.",
          text = "Please make sure the file name contains 'samples', such as samples_dataset.csv or samples.csv.",
          type = "error"
        )
        return()
      }

      # Save file
      file.copy(
        from = input$samples_csv$datapath,
        to = paste0(raw_dir(), "/samples.csv"),
        overwrite = TRUE
      )

      df <- tryCatch(
        {
          playbase::read.as_matrix(input$samples_csv$datapath)
        },
        error = function(w) {
          NULL
        }
      )
      if (is.null(df)) {
        data_error_modal(
          path = input$samples_csv$datapath,
          data_type = "samples"
        )
      } else {
        uploaded$samples.csv <- df
      }
    })

    observeEvent(input$remove_samples, {
      delete_all_files_samples <- function(value) {
        if (value) {
          uploaded$samples.csv <- NULL
          uploaded$contrasts.csv <- NULL
          checklist$samples.csv$checks <- NULL
          checklist$contrasts.csv$checks <- NULL
          checklist$samples_counts$checks <- NULL
          checklist$samples_contrasts$checks <- NULL
        }
      }

      # if contrasts is not null, warn user that it will be deleted
      if (!is.null(uploaded$contrasts.csv)) {
        shinyalert::shinyalert(
          inputId = "alert_delete_samples",
          title = "Warning",
          text = "Removing samples will also remove contrasts. Do you want to proceed?",
          type = "warning",
          showCancelButton = TRUE,
          closeOnEsc = FALSE,
          callbackR = delete_all_files_samples,
          confirmButtonText = "Yes, remove samples and contrasts",
          cancelButtonText = "Cancel"
        )
      } else {
        delete_all_files_samples(TRUE)
      }
    })

    observeEvent(input$load_example, {
      uploaded$samples.csv <- playbase::SAMPLES
    })

    observeEvent(input$create_from_samplenames, {
      counts <- uploaded$counts.csv
      samples <- playbase::createSampleInfoFromNames(colnames(counts)) 
      uploaded$samples.csv <- samples
    })

    TableModuleServer(
      "samples_datasets",
      func = table.RENDER,
      selector = "none"
    )
  })
}
