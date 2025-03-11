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

upload_table_preview_counts_server <- function(
    id,
    create_raw_dir,
    auth,
    uploaded,
    checked_matrix,
    checklist,
    scrollY,
    width,
    height,
    title,
    info.text,
    caption,
    upload_datatype) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    table_data <- shiny::reactive({
      shiny::req(!is.null(checked_matrix()))
      dt <- checked_matrix()
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


      div(
        bslib::as_fill_carrier(),
        style = "width: 100%; display: flex; ",
        if (is.null(uploaded$counts.csv)) {
          bslib::layout_columns(
            col_widths = c(-3, 6, -3),
            row_heights = list("auto", 11, 1),
            gap = "0.5rem",
            bslib::as_fill_carrier(
              bs_alert(tspan("The counts file (counts.csv) contains the gene counts for all samples. The file should be a tabular text file (.csv), where each row corresponds to a feature (i.e. genes) and each column corresponds to a sample."), closable = FALSE, translate_js = FALSE)
            ),
            bslib::card(
              if (upload_datatype() == "multi-omics") {
                ## shiny::selectInput(
                ##   ns("data_source"),
                ##   label = NULL,
                ##   choices = c("From pgx", "From csv"),
                ##   selected = "From pgx"
                ## )
                shiny::radioButtons(
                  ns("data_source"),
                  label = "Select input files from:",
                  choices = c("csv", "pgx"),
                  selected = "csv",
                  inline = TRUE
                )
              },
              if (upload_datatype() == "multi-omics") {
                shiny::conditionalPanel(
                  condition = sprintf("input['%s'] == 'pgx'", ns("data_source")),
                  div(
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
                  options = list(
                    pageLength = 5,
                    dom = "tp",
                    scrollY = TRUE
                  )
                )
              },
              if (upload_datatype() == "multi-omics") {
                shiny::conditionalPanel(
                  condition = sprintf("input['%s'] == 'csv'", ns("data_source")),
                  shiny::uiOutput(ns("dynamic_file_inputs"))#,
                )
              },
              if (upload_datatype() != "multi-omics") {
                fileInputArea(
                  ns("counts_csv"),
                  shiny::h4(tspan("Upload counts.csv", js = FALSE), class = "mb-0"),
                  multiple = FALSE,
                  accept = c(".csv"),
                  width = "100%"
                )
              },
              style = "background-color: aliceblue; border: 0.07rem dashed steelblue;"
            ),
            action_buttons
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
                bslib::nav_panel(
                  title = "Histogram",
                  br(),
                  plotOutput(ns("histogram"))
                ),
                bslib::nav_panel(
                  title = "Box plots",
                  br(),
                  plotOutput(ns("boxplots"))
                )
              )
            ),
            bslib::layout_columns(
              action_buttons,
              br(),
              uiOutput(ns("error_summary"))
            )
          )
        } ## end of if-else
      ) ## end of div
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
      DT::datatable(data.frame(
        "Dataset" = info$dataset,
        "Type" = info$datatype
      ), 
      class = "compact hover",
      rownames = FALSE,
      options = list(
        dom = 'ft',
        paging = FALSE,
        ordering = TRUE,
        info = FALSE,
        search = list(regex = FALSE, caseInsensitive = TRUE)
      ))
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

    observeEvent({
      list(input$file_input_1, input$file_input_2, input$file_input_3)
    },{
      fileinputs <- list(input$file_input_1, input$file_input_2, input$file_input_3)
      if(sum(!sapply(fileinputs,is.null))>=2) {
        shinyjs::removeClass(id = "load_selected", class = "btn-outline-primary")
        shinyjs::addClass(id = "load_selected", class = "btn-primary")
      } else {
        shinyjs::addClass(id = "load_selected", class = "btn-outline-primary")
        shinyjs::removeClass(id = "load_selected", class = "btn-primary")
      }
    })
    
    observeEvent(input$load_selected, {

      fileinputs <- list(input$file_input_1, input$file_input_2, input$file_input_3)
      numfiles <- sum(!sapply(fileinputs,is.null))
      if(numfiles < 2) {
        shinyalert::shinyalert(
          title = "Not enough input files",
          text = "Multi-omics needs minimal two data files",
          type = "error"
        )
        return(NULL)
      }
      
      col_lists <- list()
      file_names <- character()
      if (input$data_source == "pgx") { # Case from pgx not implemented yet
        #for (i in 1:length(input$available_data_table_rows_selected)) {
          info <- available_data_table()
          datasets <- info$dataset[input$available_data_table_rows_selected]
          for (i in 1:length(datasets)) {
            dataset <- datasets[i]
            pgxfile <- file.path(auth$user_dir, paste0(dataset, ".pgx"))
            df <- playbase::pgx.load(pgxfile)$counts
            col_lists[[i]] <- colnames(df)
            file_names[i] <- dataset
          }
        #}
      } else {
        for (i in 1:3) {
          file_input <- input[[paste0("file_input_", i)]]
          if (!is.null(file_input)) {
            df <- playbase::read_counts(file_input$datapath)
            col_lists[[i]] <- colnames(df)
            file_names[i] <- file_input$name
          }
        }
      }

      ## remove empty
      sel <- which(!is.na(file_names))
      col_lists <- col_lists[sel]
      file_names <- file_names[sel]
      
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
                  paste(missing, collapse=", ")
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
          col_lists[[1]] <- common_cols
        }
      }
      combined_df <- NULL
      for (i in 1:max(3, length(input$available_data_table_rows_selected))) {
        if (input$data_source == "pgx") {
          dataset <- available_data_table()[input$available_data_table_rows_selected[i], "dataset"]
          if (is.na(dataset)) {
            next
          } else {
            pgxfile <- file.path(auth$user_dir, paste0(dataset, ".pgx"))
            df <- playbase::pgx.load(pgxfile)$counts
          }
        } else {
          file_path <- input[[paste0("file_input_", i)]]$datapath
          if (is.null(file_path)) {
            next
          }
          df <- playbase::read_counts(file_path)
        }
        prefix <- switch(input[[paste0("datatype_", i)]],
          "RNA-seq" = "gx",
          "Proteomics" = "px",
          "Metabolomics" = "mx",
          "mx" # default fallback
        )
        rownames(df) <- paste0(prefix, ":", rownames(df))
        df <- df[, col_lists[[1]], drop = FALSE]
        combined_df <- rbind(combined_df, df)
      }
      uploaded$counts.csv <- combined_df
    })

    output$error_summary <- renderUI({
      div(
        style = "display: flex; justify-content: right; vertical-align: text-bottom; margin: 8px;",
        check_to_html(
          checklist$counts.csv$checks,
          pass_msg = tspan("All counts checks passed", js = FALSE),
          null_msg = tspan("Counts checks not run yet.
                            Fix any errors with counts first.", js = FALSE),
          details = FALSE
        )
      )
    })

    output$histogram <- renderPlot({
      counts <- checked_matrix()
      shiny::req(counts)
      prior <- min(counts[counts > 0], na.rm = TRUE)
      xx <- log2(prior + counts)
      # Add seed to make it deterministic
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
      counts <- checked_matrix()
      shiny::req(counts)
      prior <- min(counts[counts > 0], na.rm = TRUE)
      X <- log2(pmax(counts, 0) + prior)
      boxplot(X, ylab = tspan("counts (log2)", js = FALSE))
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
          null_msg = tspan("Counts checks not run yet.
                            Fix any errors with counts first.", js = FALSE),
          false_msg = tspan("Counts checks: warning", js = FALSE),
          details = TRUE
        )
        shinyalert::shinyalert(
          title = "Warning",
          text = err.html,
          html = TRUE
        )
      }
    })

    # pass counts to uploaded when uploaded
    observeEvent(input$counts_csv, {
      # check if counts is csv (necessary due to drag and drop of any file)
      ext <- tools::file_ext(input$counts_csv$name)
      if (!all(ext %in% c("csv", "RData"))) {
        shinyalert::shinyalert(
          title = "File format not supported.",
          text = "Please make sure the file is a CSV file.",
          type = "error"
        )
        return()
      }

      # if counts not in file name, give warning and return
      if (!any(grepl(
        "count|expression|abundance|concentration|params.rdata",
        tolower(input$counts_csv$name)
      ))) {
        shinyalert::shinyalert(
          title = tspan("Counts not in filename.", js = FALSE),
          text = tspan("Please make sure the file name contains 'counts', such as counts_dataset.csv or counts.csv.", js = FALSE),
          type = "error"
        )
        return()
      }

      # Save file
      if (!is.null(raw_dir()) && dir.exists(raw_dir())) {
        file.copy(
          from = input$counts_csv$datapath,
          to = paste0(raw_dir(), "/counts.csv"),
          overwrite = TRUE
        )
      } else { # At first raw_dir will not exist, if the user deletes and uploads a different counts it will already exist
        raw_dir(create_raw_dir(auth))
        file.copy(
          from = input$counts_csv$datapath,
          to = paste0(raw_dir(), "/counts.csv"),
          overwrite = TRUE
        )
      }

      sel <- grep("count|expression|abundance|concentration", tolower(input$counts_csv$name))
      if (length(sel)) {
        df <- tryCatch(
          {
            playbase::read_counts(input$counts_csv$datapath[sel[1]])
          },
          error = function(w) {
            NULL
          }
        )
        if (is.null(df)) {
          data_error_modal(
            path = input$counts_csv$datapath[sel[1]],
            data_type = "counts"
          )
        } else {
          uploaded$counts.csv <- df

          # if counts file contains annotation
          af <- playbase::read_annot(input$counts_csv$datapath[sel[1]])
          uploaded$annot.csv <- af
        }
      }

      sel <- grep("samples", tolower(input$counts_csv$name))
      if (length(sel)) {
        df <- tryCatch(
          {
            playbase::read_samples(input$counts_csv$datapath[sel[1]])
          },
          error = function(w) {
            NULL
          }
        )
        if (is.null(df)) {
          data_error_modal(
            path = input$counts_csv$datapath[sel[1]],
            data_type = "samples"
          )
        } else {
          uploaded$samples.csv <- df
        }
      }

      sel <- grep("contrast|comparison", tolower(input$counts_csv$name))
      if (length(sel)) {
        df <- tryCatch(
          {
            playbase::read_contrasts(input$counts_csv$datapath[sel[1]])
          },
          error = function(w) {
            NULL
          }
        )
        if (is.null(df)) {
          data_error_modal(
            path = input$counts_csv$datapath[sel[1]],
            data_type = "contrasts"
          )
        } else {
          uploaded$contrasts.csv <- df
        }
      }

      sel <- grep("params.RData", input$counts_csv$name)
      if (length(sel)) {
        if (opt$DEVMODE) {
          params <- readRDS(input$counts_csv$datapath[sel[1]])
          uploaded$samples.csv <- params$samples
          uploaded$contrasts.csv <- params$contrasts
          uploaded$counts.csv <- params$counts
        } else {
          shinyalert::shinyalert(
            title = "Error",
            text = "Invalid file: params.RData"
          )
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

    observeEvent(input$load_example, {
      uploaded$counts.csv <- playbase::COUNTS
    })


    TableModuleServer(
      "counts_datasets",
      func = table.RENDER,
      selector = "none"
    )
  }) ## end of moduleServer
} ## end of server
