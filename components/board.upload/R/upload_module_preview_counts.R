##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

## ---------------------------------------------------
## COUNTS
## ---------------------------------------------------

upload_table_preview_counts_ui <- function(id) {
  ns <- shiny::NS(id)
  uiOutput(ns("table_counts"), fill = TRUE)
}

upload_table_preview_counts_server <- function(
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
      shiny::req(!is.null(uploaded$counts.csv))
      dt <- uploaded$counts.csv
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
        DT::formatRound(columns = 1:ncol(dt), digits = 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    output$table_counts <- shiny::renderUI({
      action_buttons <- div(
        style = "display: flex; justify-content: left; margin: 8px;",
        if (is.null(uploaded$counts.csv)) {
          div(
            actionButton(
              ns("load_example"), "Load example data",
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
        ## style = "width: 100%; display: flex; justify-content: space-between; margin-bottom: 8px;",
        style = "width: 100%; display: flex; ",
        if (is.null(uploaded$counts.csv)) {
          bslib::layout_columns(
            col_widths = c(-3, 6, -3),
            row_heights = list("auto", 11, 1),
            gap = "0.5rem",
            bslib::as_fill_carrier(
              bs_alert("The counts file (counts.csv) contains the measurements (genes, proteins, etc..) for all samples. The file should be a tabular text file (.csv), where each row corresponds to a feature (i.e. genes or proteins) and each column corresponds to a sample.", closable = FALSE),
              style = "align-items: end"
            ),
            bslib::card(
              fileInputArea(
                ns("counts_csv"),
                shiny::h4("Upload abundance.csv", class = "mb-0"),
                multiple = FALSE,
                accept = c(".csv"),
                width = "100%"
              ),
              style = "background-color: aliceblue; border: 0.07rem dashed steelblue;"
            ),
            action_buttons
          )
        },
        if (!is.null(uploaded$counts.csv)) {
          bslib::layout_columns(
            col_widths = c(9, 3),
            TableModuleUI(
              ns("counts_datasets"),
              width = width,
              height = height,
              title = title,
              info.text = info.text,
              caption = caption,
              label = "",
              show.maximize = FALSE
            ),
            bslib::card(
              div(
                "Summary:",
                br(),
                check_to_html(
                  checklist$counts.csv$checks,
                  pass_msg = "All counts checks passed",
                  null_msg = "Counts checks not run yet.
                            Fix any errors with counts first."
                ),
                preview_module_legend
              )
            ),
            action_buttons
          )
        } ## end of if-else
      ) ## end of div
    })

    # pass counts to uploaded when uploaded
    observeEvent(input$counts_csv, {
      # check if counts is csv (necessary due to drag and drop of any file)
      ext <- tools::file_ext(input$counts_csv$name)[1]
      if (ext != "csv") {
        shinyalert::shinyalert(
          title = "File format not supported.",
          text = "Please make sure the file is a CSV file.",
          type = "error"
        )
        return()
      }

      # if counts not in file name, give warning and return
      if (!grepl("count", input$counts_csv$name, ignore.case = TRUE) && !grepl("expression", input$counts_csv$name, ignore.case = TRUE)) {
        shinyalert::shinyalert(
          title = "Counts not in filename.",
          text = "Please make sure the file name contains 'counts', such as counts_dataset.csv or counts.csv.",
          type = "error"
        )
        return()
      }

      df0 <- playbase::read.as_matrix(input$counts_csv$datapath)
      df <- df0
      # file.copy(df0, file.path(raw_dir(), "raw_counts.csv"))

      IS_EXPRESSION <- grepl("expression", input$counts_csv$name, ignore.case = TRUE)
      IS_COUNT <- grepl("count", input$counts_csv$name, ignore.case = TRUE)

      if (IS_EXPRESSION) {
        df <- 2**df0
        # legacy code.. it might not be always the case, maybe we should ask users if they removed zeros by adding 1 to counts... /MMM
        if (min(df0, na.rm = TRUE) > 0) df <- df - 1


        # warn user that this is happening in background
        shinyalert::shinyalert(
          title = "",
          text = "Expression counts uploaded. Converting log2 value to intensities.",
          type = "warning"
        )
      }

      uploaded$counts.csv <- df
    })

    observeEvent(input$remove_counts, {
      delete_all_files_counts <- function(value) {
        if (input$alert_delete_counts) {
          uploaded$counts.csv <- NULL
          uploaded$samples.csv <- NULL
          uploaded$contrasts.csv <- NULL
        }
      }

      # if samples is not null, warn user that it will be deleted
      if (!is.null(uploaded$samples.csv) || !is.null(uploaded$contrasts.csv)) {
        shinyalert::shinyalert(
          inputId = "alert_delete_counts",
          title = "Warning",
          text = "Removing counts will also remove samples and contrasts. Do you want to proceed?",
          type = "warning",
          showCancelButton = TRUE,
          closeOnEsc = FALSE,
          callbackR = delete_all_files_counts,
          confirmButtonText = "Remove all files",
          cancelButtonText = "Cancel"
        )

        uploaded$counts.csv <- NULL
        uploaded$samples.csv <- NULL
        uploaded$contrasts.csv <- NULL
      } else {
        uploaded$counts.csv <- NULL
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



## --------------------------------------------------------------------------
## convert list of checks to html tags for display in the data preview modal
## --------------------------------------------------------------------------

check_to_html <- function(check, pass_msg = "", null_msg = "", false_msg = "") {
  error_list <- playbase::PGX_CHECKS
  if (is.null(check)) {
    tagList(
      span(null_msg, style = "color: red"), br()
    )
  } else if (isFALSE(check)) {
    tagList(
      span(false_msg, style = "color: orange"), br()
    )
  } else {
    if (length(check) > 0) {
      tagList(
        lapply(1:length(check), function(idx) {
          error_id <- names(check)[idx]
          error_log <- check[[idx]]
          error_detail <- error_list[error_list$error == error_id, ]
          error_length <- length(error_log)
          ifelse(length(error_log) > 5, error_log <- error_log[1:5], error_log)
          if (error_detail$warning_type == "warning") {
            title_color <- "orange"
          } else if (error_detail$warning_type == "error") {
            title_color <- "red"
          }
          div(
            hr(style = "border-top: 1px solid black;"),
            span(error_detail$title, style = paste("color:", title_color)),
            br(),
            paste(error_detail$message, "\n", paste(error_length, "case(s) identified, examples:"), paste(error_log, collapse = " "), sep = " "),
            br()
          )
        }),
        hr(style = "border-top: 1px solid black;")
      )
    } else {
      if (pass_msg != "") {
        tagList(
          span(pass_msg, style = "color: green"), br()
        )
      }
    }
  }
}


preview_module_legend <- shiny::div(
  class = "pt-4",
  style = "margin-top: 150px;",
  span(style = "color: green", "Green"),
  span("= data OK. "),
  br(),
  span(style = "color: orange", "Orange"),
  span("= warning but data will still be uploaded. "),
  br(),
  span(style = "color:red", "Red"),
  span("= error and data will not be uploaded.")
)
