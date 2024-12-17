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
              bs_alert(tspan("The counts file (counts.csv) contains the gene counts for all samples. The file should be a tabular text file (.csv), where each row corresponds to a feature (i.e. genes) and each column corresponds to a sample."), closable = FALSE, translate_js = FALSE)
            ),
            bslib::card(
              fileInputArea(
                ns("counts_csv"),
                shiny::h4(tspan("Upload counts.csv", js = FALSE), class = "mb-0"),
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
              bslib::card(
                bslib::navset_pill(
                  bslib::nav_panel(
                    title = "Histogram",
                    br(),
                    plotOutput(ns("histogram"), height = "500px")
                  ),
                  bslib::nav_panel(
                    title = "Box plots",
                    br(),
                    plotOutput(ns("boxplots"), height = "500px")
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
        } ## end of if-else
      ) ## end of div
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
      xx <- log2(1 + counts)
      # Add seed to make it deterministic
      set.seed(123)
      if (nrow(xx) > 1000) xx <- xx[sample(1:nrow(xx), 1000), , drop = FALSE]
      suppressWarnings(dc <- data.table::melt(xx))
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
      X <- log2(pmax(counts, 0))
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
      if (!any(grepl("count|expression|abundance|concentration|params.rdata",
                     tolower(input$counts_csv$name)))) {
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
        df <- playbase::read_counts(input$counts_csv$datapath[sel[1]])
        uploaded$counts.csv <- df

        ## if counts file contains annotation
        af <- playbase::read_annot(input$counts_csv$datapath[sel[1]])
        uploaded$annot.csv <- af
      }

      sel <- grep("samples", tolower(input$counts_csv$name))
      if (length(sel)) {
        df <- playbase::read_samples(input$counts_csv$datapath[sel[1]])
        uploaded$samples.csv <- df
      }

      sel <- grep("contrast|comparison", tolower(input$counts_csv$name))
      if (length(sel)) {
        df <- playbase::read_contrasts(input$counts_csv$datapath[sel[1]])
        uploaded$contrasts.csv <- df
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
