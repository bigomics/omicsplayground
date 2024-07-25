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
    uploaded,
    checked_matrix,
    checklist,
    scrollY,
    width,
    height,
    title,
    info.text,
    upload_datatype,
    caption) {
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
        style = "width: 100%; display: flex; ",
        if (is.null(uploaded$counts.csv)) {
          bslib::layout_columns(
            col_widths = c(-3, 6, -3),
            row_heights = list("auto", 11, 1),
            gap = "0.5rem",
            bslib::as_fill_carrier(
              bs_alert(tspan("The counts file (counts.csv) contains the measurements (genes, proteins, etc..) for all samples. The file should be a tabular text file (.csv), where each row corresponds to a feature (i.e. genes or proteins) and each column corresponds to a sample."), closable = FALSE),
              style = "align-items: end"
            ),
            bslib::card(
              fileInputArea(
                ns("counts_csv"),
                shiny::h4(ifelse(tolower(upload_datatype()) == "proteomics", "Upload abundance.csv", "Upload counts.csv"), class = "mb-0"),
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
            col_widths = c(8, 4),
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
                br(),
                plotOutput(ns("histogram")),                       
                br(), hr(),
                "Summary:",
                br(),
                check_to_html(
                  checklist$counts.csv$checks,
                  pass_msg = tspan("All counts checks passed"),
                  null_msg = tspan("Counts checks not run yet.
                            Fix any errors with counts first."),
                  details = FALSE
                )
              ## preview_module_legend
              )
            ),
            action_buttons
          )
        } ## end of if-else
      ) ## end of div
    })

    output$histogram <- renderPlot({
      counts <- checked_matrix()
      shiny::req(counts)
      xx <- log2(1 + counts)
      if (nrow(xx) > 1000) xx <- xx[sample(1:nrow(xx), 1000), , drop = FALSE]
      suppressWarnings(dc <- data.table::melt(xx))
      dc$value[dc$value == 0] <- NA
      tt2 <- paste(nrow(counts), "genes x", ncol(counts), "samples")
      ggplot2::ggplot(dc, ggplot2::aes(x = value, color = Var2)) +
        ggplot2::geom_density() +
        ggplot2::xlab(tspan("counts (log2)")) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ggtitle(toupper(tspan("Counts")), subtitle = tt2)
    })
    
    # error pop-up alert
    observeEvent( checklist$counts.csv$checks, {
      checks <- checklist$counts.csv$checks
      if(is.null(checks)) return(NULL)
      
      if(length(checks) > 0) {
        err.html <- check_to_html(
          checks,
          pass_msg = tspan("All counts checks passed"),
          null_msg = tspan("Counts checks not run yet.
                            Fix any errors with counts first."),
          false_msg = tspan("Counts checks: warning"),
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
          title = tspan("Counts not in filename."),
          text = tspan("Please make sure the file name contains 'counts', such as counts_dataset.csv or counts.csv."),
          type = "error"
        )
        return()
      }

      df <- playbase::read_counts(input$counts_csv$datapath)

      is_expression <- grepl("expression", input$counts_csv$name, ignore.case = TRUE)
      is_count <- grepl("count", input$counts_csv$name, ignore.case = TRUE)

      if (FALSE && is_expression) {
        df <- 2**df
        if (min(df, na.rm = TRUE) >= 1) df <- df - 1

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
          text = tspan("Removing counts will also remove samples and contrasts. Do you want to proceed?"),
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


