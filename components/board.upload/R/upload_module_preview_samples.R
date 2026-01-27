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
  orig_sample_matrix,
  orig_counts_matrix,
  loaded_samples,
  sum_techreps,
  vars_selected,
  uploaded,
  checklist,
  scrollY,
  width,
  height,
  title,
  info.text,
  caption,
  upload_datatype,
  is.olink,
  public_dataset_id
) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    shiny::observe({
      if (!is.null(uploaded$counts.csv) && is.null(orig_counts_matrix())) {
        orig_counts_matrix(uploaded$counts.csv)
      }
    })
    
    shiny::observe({
      if (is.null(orig_sample_matrix()) && !is.null(uploaded$samples.csv)) {
        orig_sample_matrix(uploaded$samples.csv)
        loaded_samples(TRUE)
      }
    })

    vars_selected_pending <- reactiveVal(NULL)
    observeEvent(input$vars_selected_pending, {
      vars_selected_pending(input$vars_selected_pending)
    })

    observeEvent(input$update_vars_selected, {
      sel <- intersect(vars_selected_pending(), colnames(orig_sample_matrix()))
      if (length(sel) > 0) vars_selected(sel)
    })
    
    observeEvent(input$sum_treps, {
      shiny::req(orig_sample_matrix(), orig_counts_matrix())
      shiny::req(length(input$treps_var) > 0)
      sel <- intersect(input$treps_var, colnames(orig_sample_matrix()))
      shiny::req(length(sel) > 0)
      Y <- orig_sample_matrix()
      counts1 <- playbase::sum_treps(orig_counts_matrix(), Y[, sel[1]])
      uploaded$counts.csv <- counts1
      uploaded$samples.csv <- Y[colnames(counts1), , drop = FALSE]
      sum_techreps(TRUE)
    })

    observeEvent(input$sum_treps, {
      if (length(input$treps_var) == 0) {
        uploaded$counts.csv  <- orig_counts_matrix()
        uploaded$samples.csv <- orig_sample_matrix()
        sum_techreps(FALSE)
      }
    })

    table_data <- shiny::reactive({
      shiny::req(!is.null(uploaded$samples.csv))
      dt <- orig_sample_matrix()
      if (sum_techreps()) dt <- uploaded$samples.csv
      vars_selected <- vars_selected()
      vars_selected <- intersect(vars_selected, colnames(dt))
      dt <- dt[, vars_selected, drop = FALSE]
      uploaded$samples.csv <- dt
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
      loaded_samples(TRUE)
      return(dt)
    })

    shiny::observe({
      cols <- colnames(orig_sample_matrix())
      current_selected <- vars_selected()
      if (is.null(current_selected) || length(current_selected) == 0 || !all(current_selected %in% cols)) {
        vars_selected(cols)
      }
    })

    output$col_sel <- renderUI({
      choices <- colnames(orig_sample_matrix())
      tagList(
        checkboxGroupInput(
          ns("vars_selected_pending"),
          label = "Retain variable:",
          choices = choices,
          selected = vars_selected(),
          inline = TRUE
        ),
        actionButton(
          ns("update_vars_selected"),
          "Update",
          icon = icon("refresh"),
          class = "btn-sm btn-outline-primary"
        ),
        br()
      )
    })
    
    output$tech_rep <- renderUI({
      tagList(
        checkboxGroupInput(
          ns("treps_var"),
          label = "Technical replicates' variable:",
          choices = vars_selected(),
          selected = NULL,
          inline = TRUE
        ),
        actionButton(
          ns("sum_treps"),
          "Combine/Restore technical replicates",
          icon = icon("plus"),
          class = "btn-sm btn-outline-primary"
        ),
        br()
      )
    })

    output$add_metadata <- renderUI({
      tagList(
        withTooltip(
          checkboxInput(
            ns("add_metadata_button"),
            label = "Add metadata",
            value = FALSE
          ),
          "Expand Olink metadata by uploading an additional sample file (.csv)"
        ),
        uiOutput(ns("metadata_upload_area"))
      )
    })
    
    output$metadata_upload_area <- renderUI({
      shiny::req(input$add_metadata_button == TRUE)
      bslib::card(fileInputArea(ns("metadata_csv"),
          shiny::h4("Expand Olink metadata: upload an additional file (.csv)", class = "mb-0"),
          multiple = FALSE, accept = c(".csv"), width = "100%"),
        style = "background-color: #fffef5; border: 0.07rem dashed goldenrod;")
    })

    samples_options <- shiny::reactive({
      base_options <- shiny::tagList(
        uiOutput(ns("col_sel")),
        br(),
        uiOutput(ns("tech_rep"))
      )
      if (upload_datatype() == "proteomics" && is.olink()) {
        shiny::tagList(
          base_options,
          br(),
          uiOutput(ns("add_metadata"))
        )
      } else {
        base_options
      }
    })

    table.RENDER <- function() {
      dt <- table_data()
      req(!is.null(dt))
      DT::datatable(
        dt,
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
        div(
          if (loaded_samples()) {
            shiny::actionButton(
              ns("remove_samples"),
              "Cancel",
              icon = icon("trash-can"),
              class = "btn-sm btn-outline-danger m-1"
            )
          } else {
            shiny::actionButton(
              ns("load_example"), "Load example",
              class = "btn-sm btn-outline-primary m-1"
            )
          }
        ),
        div(
          shiny::actionButton(
            ns("check_documentation_samples"),
            "Read documentation",
            class = "btn-sm btn-outline-primary m-1",
            onclick = "window.open('https://omicsplayground.readthedocs.io/en/latest/dataprep/samples/', '_blank')"
          )
        )
      )

      div(
        bslib::as_fill_carrier(),
        if (!loaded_samples() && public_dataset_id() == "") {
          bslib::layout_columns(
            col_widths = c(-3, 6, -3),
            row_heights = list("auto", 8, 1, 2),
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
            action_buttons,
            br()
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
                options = samples_options(),
                label = "",
                show.maximize = FALSE
              ),
              bslib::navset_card_pill(
                bslib::nav_panel(
                  title = "UMAP",
                  br(),
                  plotOutput(ns("umap"))
                ),
                bslib::nav_panel(
                  title = "Distribution",
                  br(),
                  plotOutput(ns("phenotype_stats"))
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
          pass_msg = "Samples matrix OK. ",
          null_msg = "Samples checks not run yet. Fix any errors with samples first. ",
          false_msg = "Samples checks: warning ",
          details = FALSE
        ),
        check_to_html(
          checklist$samples_counts$checks,
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
      shiny::req(nrow(Y))
      shiny::validate(shiny::need(ncol(Y) > 0, "Please select at least 1 variable."))
      sel <- grep("group|condition", colnames(Y), ignore.case = TRUE)
      sel <- head(c(sel, 1), 1)
      y <- Y[, sel]
      # Fix: handle case where sel col is NA (plot errors)
      if (all(is.na(y))) {
        non_na_cols <- which(colSums(!is.na(Y)) > 0)
        if (length(non_na_cols) > 0) {
          y <- Y[, non_na_cols[1]]
        }
      }
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
          shinyalert::shinyalert(title = "Warning", text = err.html, html = TRUE)
        }
      }
    )

    ## pass counts to uploaded when uploaded
    observeEvent(input$samples_csv, {
      counts <- uploaded$counts.csv
      shiny::req(nrow(counts))
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
      datafile <- input$samples_csv$datapath
      file.copy(from = datafile, to = paste0(raw_dir(), "/samples.csv"), overwrite = TRUE)
      df <- tryCatch({ playbase::read.as_matrix(datafile) },
        error = function(w) { NULL })      
      if (is.null(df)) {
        data_error_modal(path = datafile, data_type = "samples")
      } else {
        uploaded$samples.csv <- df
      }
      orig_sample_matrix(df)
      loaded_samples(TRUE)
    })

    observeEvent(input$metadata_csv, {
      shiny::req(input$metadata_csv)
      c1 <- (tools::file_ext(input$metadata_csv$name)[1] != "csv")
      c2 <- (!grepl("sample", input$metadata_csv$name, ignore.case = TRUE))
      if (c1 | c2) {
        shinyalert::shinyalert(
          title = "File format not supported.",
          text = "Please make sure the file is a CSV file and contains 'samples' in the file name.",
          type = "error"
        )
        return()
      }
      datafile <- input$metadata_csv$datapath
      new_samples <- tryCatch({ playbase::read.as_matrix(datafile) },
        error = function(w) { NULL })
      if (is.null(new_samples)) {
        data_error_modal(path = datafile, data_type = "metadata")
        return()
      }
      samples <- uploaded$samples.csv
      shiny::req(samples)
      cm <- intersect(as.character(rownames(samples)), as.character(rownames(new_samples)))
      if (length(cm) == 0) {
        shinyalert::shinyalert(
          title = "No matching samples.",
          text = "The newly uploaded metadata file does not share any sample identifiers with the current samples",
          type = "error"
        )
        return()
      }
      if ((length(cm) != nrow(samples)) | (length(cm) != nrow(new_samples))) {
        shinyalert::shinyalert(
          title = "Samples mismatch.",
          text = "The newly uploaded sample file contains a different set of samples than the current metadata. We will intersect these.",
          type = "error"
        )
        samples <- samples[cm, ]
        new_samples <- new_samples[cm, ]
      }
      new_cols <- setdiff(colnames(new_samples), colnames(samples))
      if (length(new_cols) > 0) {
        merged <- cbind(samples, new_samples[rownames(samples), new_cols, drop = FALSE])
        uploaded$samples.csv <- merged
        orig_sample_matrix(merged)
        vars_selected(colnames(merged))
        shinyalert::shinyalert(
          title = "Metadata added",
          text = paste("Added", length(new_cols), "new variables:", paste(new_cols, collapse = ", ")),
          type = "success"
        )
      } else {
        shinyalert::shinyalert(
          title = "No new columns added",
          text = "The newly upload metadata file does not contain any new variables.",
          type = "warning"
        )
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
      loaded_samples(FALSE)
      vars_selected(NULL)
      orig_sample_matrix(NULL)
    })

    observeEvent(input$load_example, {
      df <- playbase::SAMPLES
      if (upload_datatype() == "multi-omics") df <- playbase::SAMPLES_MO
      if (upload_datatype() == "scRNA-seq") df <- playdata::GSE243639_scRNAseq_samples
      uploaded$samples.csv <- df
      loaded_samples(TRUE)
      orig_sample_matrix(df)
    })

    TableModuleServer(
      "samples_datasets",
      func = table.RENDER,
      selector = "none"
    )
  })
}
