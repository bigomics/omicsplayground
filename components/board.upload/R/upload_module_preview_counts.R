##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

## ---------------------------------------------------
## COUNTS UPLOAD (for wizard dialog)
## ---------------------------------------------------

upload_table_preview_counts_ui <- function(id) {
  ns <- shiny::NS(id)
  uiOutput(ns("table_counts"), fill = TRUE)
}

upload_table_preview_counts_server <- function(id,
                                               create_raw_dir,
                                               auth,
                                               uploaded,
                                               checked_matrix,
                                               is_logscale,
                                               checklist,
                                               scrollY,
                                               width,
                                               height,
                                               title,
                                               info.text,
                                               caption,
                                               upload_datatype,
                                               is.olink,
                                               public_dataset_id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    GEO_alert_shown <- reactiveVal(FALSE)

    table_data <- shiny::reactive({
      shiny::req(!is.null(uploaded$counts.csv))
      dt <- uploaded$counts.csv
      nrow0 <- nrow(dt)
      ncol0 <- ncol(dt)
      MAXROW <- 1000
      MAXCOL <- 20
      if (nrow(dt) > MAXROW) {
        dt <- head(dt, MAXROW)
        dt <- rbind(dt, rep(NA, ncol(dt)))
        n1 <- nrow0 - MAXROW
        rownames(dt)[nrow(dt)] <- paste0("[+", n1, " rows]")
      }
      if (ncol(dt) > MAXCOL) {
        dt <- dt[, 1:MAXCOL]
        dt <- cbind(dt, rep(NA, nrow(dt)))
        n1 <- ncol0 - MAXCOL
        colnames(dt)[ncol(dt)] <- paste0("[+", n1, " columns]")
      }
      dt
    })

    table.RENDER <- function() {
      dt <- table_data()
      req(!is.null(dt))
      is.integer <- is.integer(dt) || all(round(dt) == dt, na.rm = TRUE)
      digits <- ifelse(is.integer, 0, 2)
      DT::datatable(dt,
        class = "compact hover",
        rownames = TRUE,
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = 1),
        options = list(
          dom = "lfrtip",
          scroller = TRUE,
          scrollX = TRUE,
          scrollY = scrollY,
          scrollResize = TRUE,
          deferRender = TRUE
        )
      ) %>%
        DT::formatRound(columns = 1:ncol(dt), digits = digits) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    output$table_counts <- shiny::renderUI({
      action_buttons <- div(
        style = "display: flex; justify-content: left; margin: 8px;",
        if (is.null(uploaded$counts.csv)) {
          div(
            if (upload_datatype() == "multi-omics") {
              actionButton(
                ns("load_selected"), "Upload selected files",
                class = "btn-sm btn-outline-primary m-1"
              )
            },
            actionButton(
              ns("load_example"), "Load example",
              class = "btn-sm btn-outline-primary m-1"
            )
          )
        } else {
          div(
            shiny::actionButton(
              ns("remove_counts"),
              "Cancel",
              icon = icon("trash-can"),
              class = "btn-sm btn-outline-danger m-1"
            )
          )
        },
        shiny::actionButton(
          ns("check_documentation_counts"),
          "Read documentation",
          class = "btn-sm btn-outline-primary m-1",
          onclick = "window.open('https://omicsplayground.readthedocs.io/en/latest/dataprep/counts/', '_blank')"
        )
      )

      if (public_dataset_id() != "") {
        ID <- public_dataset_id()
        msg <- paste0("Retrieving ", ID, " from GEO, ReCount, or ArrayExpress...<br>", "Please wait. Most datasets take 2-3 mins.")

        showModal(modalDialog(
          div(
            id = "custom-progress-modal", HTML(msg),
            div(
              id = "custom-progress-container",
              div(id = "custom-progress-bar")
            )
          ),
          footer = NULL, fade = FALSE
        ))
        GEO <- tryCatch(
          {
            playbase::pgx.getGEOseries(accession = ID, archs.h5 = NULL, get.info = FALSE)
          },
          error = function(w) {
            NULL
          }
        )
        removeModal()

        if (!is.null(GEO)) {
          if (!GEO_alert_shown()) {
            msg <- paste0("Success! ", ID, " found in ", GEO[["source"]], ".\nWe're preparing it...")
            GEO_alert_shown(TRUE)
          }
          uploaded$counts.csv <- GEO[["counts"]]
          uploaded$samples.csv <- GEO[["samples"]]
          cm <- intersect(colnames(GEO[["counts"]]), rownames(GEO[["samples"]]))
          GEO <- NULL
        }

        if (!is.null(uploaded$counts.csv)) {
          bslib::layout_columns(
            col_widths = 12,
            bslib::layout_columns(
              col_widths = c(8, 4),
              TableModuleUI(
                ns("counts_datasets"),
                width = width,
                height = height,
                title = title,
                info.text = info.text,
                caption = caption,
                label = "",
                show.maximize = FALSE,
                translate_js = FALSE
              ),
              bslib::navset_card_pill(
                bslib::nav_panel(title = "Histogram", br(), plotOutput(ns("histogram"))),
                bslib::nav_panel(title = "Box plots", br(), plotOutput(ns("boxplots")))
              )
            ),
            bslib::layout_columns(action_buttons, br(), uiOutput(ns("error_summary")))
          )
        } else {
          div("No counts data available for this public dataset.")
        }
      } else {
        div(
          bslib::as_fill_carrier(),
          style = "width: 100%; display: flex; ",
          if (is.null(uploaded$counts.csv)) {
            if (upload_datatype() == "proteomics") {
              msg <- "The counts file (counts.csv) contains the gene counts for all samples. For proteomics data types other than Olink NPX, the file should be a tabular text file (.csv), where each row corresponds to a feature (i.e. genes) and each column corresponds to a sample. For Olink NPX, the uploaded file needs to be the standard Olink format and can be a parquet file."
            } else {
              msg <- "The counts file (counts.csv) contains the gene counts for all samples. The file should be a tabular text file (.csv), where each row corresponds to a feature (i.e. genes) and each column corresponds to a sample."
            }
            bslib::layout_columns(
              col_widths = c(-3, 6, -3),
              row_heights = list("auto", 8, 1, 2),
              gap = "0.5rem",
              bslib::as_fill_carrier(bs_alert(tspan(msg), closable = FALSE, translate_js = FALSE)),
              bslib::card(
                if (upload_datatype() == "multi-omics") {
                  shiny::radioButtons(
                    ns("data_source"),
                    label = "Select input files from:",
                    choices = c("multi-csv", "pgx", "single-csv"),
                    selected = "multi-csv",
                    inline = TRUE
                  )
                },
                if (upload_datatype() == "multi-omics") {
                  shiny::conditionalPanel(
                    condition = sprintf("input['%s'] == 'pgx'", ns("data_source")),
                    div(
                      div(
                        style = "margin-bottom: 15px;",
                        shiny::radioButtons(
                          ns("data_processing"),
                          label = "Select data processing level:",
                          choices = c("Raw" = "raw", "Normalized (also batch corrected if selected on computation)" = "normalized"),
                          selected = "raw",
                          inline = TRUE
                        )
                      ),
                      div(
                        style = "display: flex; align-items: center; gap: 10px; margin-bottom: 10px;",
                        span("Selected:", style = "font-weight: bold;"),
                        textOutput(ns("selected_rows_text")),
                      ),
                      div(
                        style = "height: 350px; overflow-y: auto;",
                        DT::DTOutput(ns("available_data_table"))
                      ),
                      style = "width: 100%; max-height: 400px; overflow-y: auto;"
                    ),
                    selection = "multiple",
                    options = list(pageLength = 5, dom = "tp", scrollY = TRUE)
                  )
                },
                if (upload_datatype() == "multi-omics") {
                  shiny::conditionalPanel(
                    condition = sprintf("input['%s'] == 'multi-csv'", ns("data_source")),
                    shiny::uiOutput(ns("dynamic_file_inputs")) # ,
                  )
                },
                if (upload_datatype() == "multi-omics") {
                  shiny::conditionalPanel(
                    condition = sprintf("input['%s'] == 'single-csv'", ns("data_source")),
                    fileInputArea(
                      ns("counts_csv"),
                      shiny::h4(tspan("Upload counts.csv", js = FALSE), class = "mb-0"),
                      multiple = FALSE,
                      accept = c(".csv"),
                      width = "100%"
                    )
                  )
                },
                if (upload_datatype() != "multi-omics") {
                  if (upload_datatype() == "scRNA-seq") {
                    fileInputArea(
                      ns("counts_csv"),
                      shiny::h4(tspan("Upload counts.csv/.h5", js = FALSE), class = "mb-0"),
                      multiple = FALSE,
                      accept = c(".csv", ".h5"),
                      width = "100%"
                    )
                  } else {
                    fileInputArea(
                      ns("counts_csv"),
                      shiny::h4(tspan("Upload counts.csv", js = FALSE), class = "mb-0"),
                      multiple = FALSE,
                      accept = if (upload_datatype() == "proteomics" && is.olink()) c(".csv", ".parquet") else c(".csv"),
                      width = "100%"
                    )
                  }
                },
                style = "background-color: aliceblue; border: 0.07rem dashed steelblue;"
              ),
              action_buttons,
              br()
            )
          },
          if (!is.null(uploaded$counts.csv)) {
            bslib::layout_columns(
              col_widths = 12,
              bslib::layout_columns(
                col_widths = c(8, 4),
                TableModuleUI(
                  ns("counts_datasets"),
                  width = width,
                  height = height,
                  title = title,
                  info.text = info.text,
                  caption = caption,
                  label = "",
                  show.maximize = FALSE,
                  translate_js = FALSE
                ),
                bslib::navset_card_pill(
                  bslib::nav_panel(title = "Histogram", br(), plotOutput(ns("histogram"))),
                  bslib::nav_panel(title = "Box plots", br(), plotOutput(ns("boxplots")))
                )
              ),
              bslib::layout_columns(action_buttons, br(), uiOutput(ns("error_summary")))
            )
          } ## end of if-else
        ) ## end of div
      }
    })

    output$selected_rows_text <- renderText({
      info <- available_data_table()
      paste(info$dataset[input$available_data_table_rows_selected], collapse = ", ")
    })

    available_data_table <- reactive({
      pgxdir <- auth$user_dir
      info <- playbase::pgxinfo.read(pgxdir, file = "datasets-info.csv")
      info <- info[info$datatype != "multi-omics", ]
      info
    })

    output$available_data_table <- DT::renderDT({
      info <- available_data_table()
      DT::datatable(data.frame("Dataset" = info$dataset, "Type" = info$datatype),
        class = "compact hover",
        rownames = FALSE,
        options = list(
          dom = "ft",
          paging = FALSE,
          ordering = TRUE,
          info = FALSE,
          search = list(regex = FALSE, caseInsensitive = TRUE)
        )
      )
    })

    output$dynamic_file_inputs <- renderUI({
      bslib::layout_column_wrap(
        style = bslib::css(grid_template_columns = "8fr 3fr 1fr"),
        class = "m-0",
        fileInput(ns("file_input_1"), label = NULL, multiple = FALSE, accept = c(".csv")),
        selectInput(ns("datatype_1"), label = NULL, choices = c("RNA-seq", "Proteomics", "Metabolomics")),
        actionButton(ns("remove_input_1"), label = NULL, icon = icon("xmark"), class = "btn-sm btn-outline-danger"),
        fileInput(ns("file_input_2"), label = NULL, multiple = FALSE, accept = c(".csv")),
        selectInput(ns("datatype_2"), label = NULL, choices = c("RNA-seq", "Proteomics", "Metabolomics"), selected = "Proteomics"),
        actionButton(ns("remove_input_2"), label = NULL, icon = icon("xmark"), class = "btn-sm btn-outline-danger"),
        fileInput(ns("file_input_3"), label = NULL, multiple = FALSE, accept = c(".csv")),
        selectInput(ns("datatype_3"), label = NULL, choices = c("RNA-seq", "Proteomics", "Metabolomics"), selected = "Metabolomics"),
        actionButton(ns("remove_input_3"), label = NULL, icon = icon("xmark"), class = "btn-sm btn-outline-danger")
      )
    })

    lapply(1:3, function(i) {
      observeEvent(input[[paste0("remove_input_", i)]], {
        shinyjs::hide(paste0("file_input_", i))
        shinyjs::hide(paste0("datatype_", i))
        shinyjs::hide(paste0("remove_input_", i))
      })
    })

    observeEvent(
      {
        list(input$file_input_1, input$file_input_2, input$file_input_3)
      },
      {
        fileinputs <- list(input$file_input_1, input$file_input_2, input$file_input_3)
        if (sum(!sapply(fileinputs, is.null)) >= 2) {
          shinyjs::removeClass(id = "load_selected", class = "btn-outline-primary")
          shinyjs::addClass(id = "load_selected", class = "btn-primary")
        } else {
          shinyjs::addClass(id = "load_selected", class = "btn-outline-primary")
          shinyjs::removeClass(id = "load_selected", class = "btn-primary")
        }
      }
    )

    observeEvent(input$load_selected, {
      # Validate minimum number of inputs based on data source
      enough_inputs <- TRUE
      error_title <- ""

      if (input$data_source == "pgx") {
        if (length(input$available_data_table_rows_selected) < 2) {
          enough_inputs <- FALSE
          error_title <- "Not enough datasets selected"
        }
      } else if (input$data_source == "multi-csv") {
        fileinputs <- list(input$file_input_1, input$file_input_2, input$file_input_3)
        numfiles <- sum(!sapply(fileinputs, is.null))
        if (numfiles < 2) {
          enough_inputs <- FALSE
          error_title <- "Not enough input files"
        }
      }

      if (!enough_inputs) {
        shinyalert::shinyalert(
          title = error_title,
          text = "Multi-omics needs minimal two data sources",
          type = "error"
        )
        return(NULL)
      }

      # Load data once and cache both data frames and column names
      samples_cache <- list()
      data_cache <- list()
      col_lists <- list()
      file_names <- character()
      datatypes <- character()

      if (input$data_source == "pgx") {
        info <- available_data_table()
        datasets <- info$dataset[input$available_data_table_rows_selected]
        for (i in 1:length(datasets)) {
          dataset <- datasets[i]
          pgxfile <- file.path(auth$user_dir, paste0(dataset, ".pgx"))
          pgx_data <- playbase::pgx.load(pgxfile)
          # Use raw or normalized data based on selection
          if (input$data_processing == "normalized" && !is.null(pgx_data$X)) {
            df <- pgx_data$X
          } else {
            df <- pgx_data$counts
          }
          samples_cache[[i]] <- pgx_data$samples
          data_cache[[i]] <- df
          col_lists[[i]] <- colnames(df)
          file_names[i] <- dataset
          datatypes[i] <- info$datatype[input$available_data_table_rows_selected[i]]
        }
      } else {
        for (i in 1:3) {
          file_input <- input[[paste0("file_input_", i)]]
          if (!is.null(file_input)) {
            df <- playbase::read_counts(file_input$datapath)
            data_cache[[i]] <- df
            col_lists[[i]] <- colnames(df)
            file_names[i] <- file_input$name
            datatypes[i] <- input[[paste0("datatype_", i)]]
          }
        }
      }

      ## remove empty
      sel <- which(!is.na(file_names))
      data_cache <- data_cache[sel]
      col_lists <- col_lists[sel]
      file_names <- file_names[sel]
      datatypes <- datatypes[sel]

      if (length(col_lists) > 1) {
        common_cols <- Reduce(intersect, col_lists)
        if (length(common_cols) == 0) {
          shinyalert::shinyalert(
            title = "No Common Columns",
            text = "The uploaded files have no columns in common",
            type = "error"
          )
          return(NULL)
        } else {
          all_cols <- unique(unlist(col_lists))
          excluded_cols <- setdiff(all_cols, common_cols)
          if (length(excluded_cols) > 0) {
            mismatch_text <- NULL
            for (i in seq_along(col_lists)) {
              missing <- setdiff(all_cols, col_lists[[i]])
              if (length(missing) > 0) {
                mismatch_text <- paste0(
                  mismatch_text,
                  "\nFile '", file_names[i], "' is missing columns: ",
                  paste(missing, collapse = ", ")
                )
              }
            }
            shinyalert::shinyalert(
              title = "Using Common Columns Only",
              text = paste0(
                "Some columns were not present in all files and will be excluded.",
                mismatch_text,
                "\n\nProceeding with ", length(common_cols), " common columns."
              ),
              type = "warning"
            )
          }
        }
      } else {
        common_cols <- col_lists[[1]]
      }

      # Use cached data instead of reloading
      combined_df <- NULL
      for (i in 1:length(data_cache)) {
        df <- data_cache[[i]]
        dt <- tolower(datatypes[i])
        prefix <- "gx"
        if (grepl("proteomics", dt)) prefix <- "px"
        if (grepl("metabolomics|lipidomics", dt)) prefix <- "mx"
        if (grepl("microarray|micro.array|rna|rnatranscriptomics", dt)) prefix <- "gx"
        if (grepl("mirna|mi.rna", dt)) prefix <- "mi"
        rownames(df) <- paste0(prefix, ":", rownames(df))
        df <- df[, common_cols, drop = FALSE]
        combined_df <- rbind(combined_df, df)
      }

      # Combine samples_cache by intersecting rownames (only for pgx multi-omics upload)
      if (length(samples_cache) > 0) {
        # Find common rownames in all sample files (no prefixing!)
        rownames_list <- lapply(samples_cache, rownames)
        common_rows <- Reduce(intersect, rownames_list)
        combined_samples <- samples_cache[[1]][common_rows, , drop = FALSE]
        uploaded$samples.csv <- combined_samples
      }

      uploaded$counts.csv <- combined_df
    })

    output$error_summary <- renderUI({
      div(
        style = "display: flex; justify-content: right; vertical-align: text-bottom; margin: 8px;",
        check_to_html(
          checklist$counts.csv$checks,
          pass_msg = tspan("All counts checks passed", js = FALSE),
          null_msg = tspan("Counts checks not run yet. Fix any errors with counts first.", js = FALSE),
          details = FALSE
        )
      )
    })

    output$histogram <- renderPlot({
      counts <- uploaded$counts.csv
      shiny::req(counts)
      xx <- counts
      if (!is_logscale()) {
        prior <- min(counts[counts > 0], na.rm = TRUE)
        xx <- log2(prior + counts)
      }
      set.seed(123)
      if (nrow(xx) > 1000) xx <- xx[sample(1:nrow(xx), 1000), , drop = FALSE]
      suppressWarnings(dc <- reshape2::melt(xx))
      dc$value[dc$value == 0] <- NA
      tt2 <- paste(nrow(counts), tspan("genes x", js = FALSE), ncol(counts), "samples")
      ggplot2::ggplot(dc, ggplot2::aes(x = value, color = Var2)) +
        ggplot2::geom_density() +
        ggplot2::xlab(tspan("counts (log2)", js = FALSE)) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ggtitle(toupper(tspan("Counts", js = FALSE)), subtitle = tt2)
    })

    output$boxplots <- renderPlot({
      counts <- uploaded$counts.csv
      shiny::req(counts)
      xx <- counts
      if (!is_logscale()) {
        prior <- min(xx[xx > 0], na.rm = TRUE)
        xx <- log2(pmax(xx, 0) + prior)
      }
      # Downsample to 40 columns as we do on qc/bc tab
      if (ncol(xx) > 40) xx <- xx[, sample(1:ncol(xx), 40)]
      boxplot(xx, ylab = tspan("counts (log2)", js = FALSE))
    })

    # error pop-up alert
    observeEvent(checklist$counts.csv$checks, {
      checks <- checklist$counts.csv$checks
      if (is.null(checks)) {
        return(NULL)
      }
      if (length(checks) > 0) {
        err.html <- check_to_html(
          checks,
          pass_msg = tspan("All counts checks passed", js = FALSE),
          null_msg = tspan("Counts checks not run yet. Fix any errors with counts first.", js = FALSE),
          false_msg = tspan("Counts checks: warning", js = FALSE),
          details = TRUE
        )
        shinyalert::shinyalert(title = "Warning", text = err.html, html = TRUE)
      }
    })

    # pass counts to uploaded when uploaded
    observeEvent(input$counts_csv, {
      ext <- tools::file_ext(input$counts_csv$name)
      dtypes <- c("RNA-seq", "mRNA microarray", "proteomics", "metabolomics", "lipidomics")
      c1 <- (!(upload_datatype() %in% dtypes && ext %in% c("csv", "RData")))
      c2 <- (!(upload_datatype() == "scRNA-seq" && ext %in% c("csv", "h5")))
      c3 <- (!(upload_datatype() == "proteomics" && is.olink() && ext %in% c("csv", "parquet"))) 
      if (c1 & c2 & c3) {
        shinyalert::shinyalert(
          title = "File format not supported.",
          text = "Please upload a .csv file. For scRNA-seq, h5 format is allowed. For Olink NPX data, parquet format is allowed.",
          type = "error"
        )
        return()
      }

      # if counts not in file name, give warning and return
      ss <- "count|expression|abundance|concentration|params.rdata"
      if (!any(grepl(ss, tolower(input$counts_csv$name)))) {
        shinyalert::shinyalert(
          title = tspan("Counts not in filename.", js = FALSE),
          text = tspan("Please ensure the file name contains 'counts', e.g., counts_dataset.csv or counts.csv.", js = FALSE),
          type = "error"
        )
        return()
      }

      # Save file
      # At first raw_dir will not exist, if the user deletes and uploads a different counts it will already exist
      if (!is.null(raw_dir()) && dir.exists(raw_dir())) {
        file.copy(
          from = input$counts_csv$datapath,
          to = paste0(raw_dir(), "/counts.csv"),
          overwrite = TRUE
        )
      } else {
        raw_dir(create_raw_dir(auth))
        file.copy(
          from = input$counts_csv$datapath,
          to = paste0(raw_dir(), "/counts.csv"),
          overwrite = TRUE
        )
      }

      ## ---counts---##
      sel <- grep("count|expression|abundance|concentration", tolower(input$counts_csv$name))
      if (length(sel)) {

        datafile <- input$counts_csv$datapath[sel[1]]
        datafile.name <- input$counts_csv$name
        file.ext <- tools::file_ext(datafile.name)

        if (upload_datatype() == "scRNA-seq" && file.ext == "h5") {
          df <- tryCatch({ playbase::read_h5_counts(datafile) }, error = function(w) { NULL } )
          if (is.null(df)) {
            shinyalert::shinyalert(
              title = "Error",
              text = "Error: there may be an issue with the uploaded h5 format. Please fix it & re-upload.",
              type = "error"
            )
          }
        } else {
          df.samples <- NULL
          if (upload_datatype() == "proteomics" && is.olink()) {
            df0 <- tryCatch(
            {
              playbase::read_Olink_NPX(datafile)
            },
            error = function(w) {
              NULL
            }
            )
            if (is.null(df0)) {
              shinyalert::shinyalert(
                title = "Error",
                text = "Your data may not be in Official Olink format. Please check."
              )
            } else {              
              df <- df0[["counts"]]
              df.samples <- df0[["samples"]]
              rm(df0)
            }
          } else if (upload_datatype() == "proteomics" && !is.olink()) {
            df <- tryCatch({ playbase::read_counts(datafile) }, error = function(w) { NULL } )
            if (is.null(df)) {
              df <- tryCatch({ playbase::read_spectronaut(datafile) }, error = function(w) { NULL } )
              if (!is.null(df)) {
                char.cols <- which(sapply(df, class) == "character")
                if (length(char.cols) > 0) {
                  uploaded$annot.csv <- df[, names(char.cols), drop = FALSE]
                  df <- df[, colnames(df) != names(char.cols), drop = FALSE]
                  df <- as.matrix(df)
                }
              }
            }
          } else {
            df <- tryCatch({ playbase::read_counts(datafile) }, error = function(w) { NULL } )
          }
        }
      } else {
        df <- tryCatch(
        {
          playbase::read_counts(datafile)
        },
        error = function(w) {
          NULL
        }
        )
      }
      
      
      file.ext <- tools::file_ext(input$counts_csv$name)
      if (is.null(df) & file.ext != "h5") {
        data_error_modal(path = datafile, data_type = "counts")
      } else {
        uploaded$counts.csv <- df
        if (is.null(uploaded$annot.csv)) {
          ann.data <- if (! file.ext %in% c("h5", "parquet")) playbase::read_annot(datafile) else NULL
          uploaded$annot.csv <- ann.data
        }
      }

      if (upload_datatype() == "proteomics" && is.olink() && !is.null(df.samples)) {
        uploaded$samples.csv <- df.samples
      }

      sel <- grep("params.RData", input$counts_csv$name)
      if (length(sel)) {
        if (opt$DEVMODE) {
          params <- readRDS(input$counts_csv$datapath[sel[1]])
          uploaded$samples.csv <- params$samples
          uploaded$contrasts.csv <- params$contrasts
          uploaded$counts.csv <- params$counts
        } else {
          shinyalert::shinyalert(title = "Error", text = "Invalid file: params.RData")
        }
      }
    }) ## end observeEvent input$counts.csv


    observeEvent(input$remove_counts, {
      delete_all_files_counts <- function(value) {
        if (value) {
          uploaded$counts.csv <- NULL
          uploaded$samples.csv <- NULL
          uploaded$contrasts.csv <- NULL
          checklist$counts.csv$checks <- NULL
          checklist$samples.csv$checks <- NULL
          checklist$contrasts.csv$checks <- NULL
          checklist$samples_counts$checks <- NULL
          checklist$samples_contrasts$checks <- NULL
        }
      }

      # if samples is not null, warn user that it will be deleted
      if (!is.null(uploaded$samples.csv) || !is.null(uploaded$contrasts.csv)) {
        shinyalert::shinyalert(
          inputId = "alert_delete_counts",
          title = "Warning",
          text = tspan("Removing counts will also remove samples and contrasts. Do you want to proceed?", js = FALSE),
          type = "warning",
          showCancelButton = TRUE,
          closeOnEsc = FALSE,
          callbackR = delete_all_files_counts,
          confirmButtonText = "Remove all files",
          cancelButtonText = "Cancel"
        )
      } else {
        delete_all_files_counts(TRUE)
      }
    })

    shiny::observeEvent(input$load_example, {
      if (upload_datatype() == "multi-omics") {
        uploaded$counts.csv <- playbase::COUNTS_MO
      } else if (is.olink()) {
        shinyalert::shinyalert(
          title = "Error",
          text = "Olink NPX example data not yet available. Please upload yours or change data type."
        )
      } else if (upload_datatype() == "scRNA-seq") {
        uploaded$counts.csv <- playdata::GSE243639_scRNAseq_counts
      } else {
        uploaded$counts.csv <- playbase::COUNTS
      }
    })

    TableModuleServer(
      "counts_datasets",
      func = table.RENDER,
      selector = "none"
    )
  }) ## end of moduleServer
} ## end of server
