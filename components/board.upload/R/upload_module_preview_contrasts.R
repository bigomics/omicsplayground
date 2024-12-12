##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##


upload_table_preview_contrasts_ui <- function(id) {
  ns <- shiny::NS(id)
  uiOutput(ns("table_contrasts"), fill = TRUE)
}

upload_table_preview_contrasts_server <- function(
    id,
    uploaded,
    checklist,
    scrollY,
    width,
    height,
    title,
    info.text,
    caption,
    checked_counts,
    checked_samples,
    checked_contrasts,
    show_comparison_builder,
    selected_contrast_input,
    upload_wizard) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns


    table_data <- shiny::reactive({
      shiny::req(!is.null(uploaded$contrasts.csv))
      dt <- uploaded$contrasts.csv
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
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    output$table_contrasts <- shiny::renderUI({
      # I will set selected_contrast_input to true for now (selection page will be disabled)
      # if(selected_contrast_input() == FALSE) {
      # ask user if preferrence is upload contrast or create contrast online
      # TODO remove this is not needed
      # div(
      #   style = "display: flex; gap: 80px; flex-direction: column; justify-content: center; align-items: center; height: 50%;",
      #   div(h4("Please choose one of the following options:")),
      #   div(
      #     style = "display: flex; justify-content: center; gap: 20px;",
      #     actionButton(
      #     ns("goUploadComparison"),
      #     label = " Upload comparisons file",
      #     class = "btn btn-primary m-1",
      #     icon("upload")
      #     ),
      #   actionButton(
      #     ns("goOnlineComparison"),
      #     label = " Build my comparisons",
      #     class = "btn btn-primary m-1",
      #     icon("pen-to-square")
      #     )
      #   )

      # )
      # } else {
      # only display buttons if goOnlineComparison is false

      action_buttons1 <- div(
        style = "display: flex; justify-content: left; margin-bottom: 20px;",
        actionButton(
          ns("goUploadComparison"),
          label = "Upload comparisons file",
          class = "btn-sm btn-secondary m-1",
          icon("upload")
        ),
        withTooltip(
          shiny::actionButton(
            ns("autocontrast"),
            "Auto-detect comparisons",
            ## icon = icon("plus"),
            class = "btn-sm btn-outline-primary  m-1"
          ),
          "If you are feeling lucky, try this to automatically create comparisons.",
          placement = "top", options = list(container = "body")
        ),
        actionButton(
          ns("check_documentation"),
          "Read documentation",
          class = "btn-sm btn-outline-primary  m-1",
          onclick = "window.open('https://omicsplayground.readthedocs.io/en/latest/dataprep/contrasts/', '_blank')"
        )
      )

      action_buttons2 <- div(
        style = "display: flex; justify-content: left; margin-bottom: 8px;",
        if (!is.null(uploaded$contrasts.csv)) {
          div(shiny::actionButton(
            ns("remove_contrasts"),
            "Cancel",
            icon = icon("trash-can"),
            class = "btn-sm btn-outline-danger m-1"
          ))
        } else {
          div(
            actionButton(
              ns("run_build_comparisons"), "Build my comparisons",
              class = "btn-sm btn-secondary m-1",
              icon("pen-to-square")
            ),
            actionButton(
              ns("load_example"), "Load example data",
              class = "btn-sm btn-outline-primary m-1"
            )
          )
        },
        actionButton(
          ns("check_documentation_contrasts"),
          "Read documentation",
          class = "btn-sm btn-outline-primary m-1",
          onclick = "window.open('https://omicsplayground.readthedocs.io/en/latest/dataprep/contrasts/', '_blank')"
        )
      )

      comparison_builder_ui <-
        div(
          bslib::as_fill_carrier(),
          if (is.null(uploaded$contrasts.csv)) {
            bslib::layout_columns(
              col_widths = c(-3, 6, -3),
              row_heights = list("auto", 11, 1),
              gap = "0.3rem",
              bslib::as_fill_carrier(
                bs_alert("The comparison file (comparisons.csv) is an optional input file. The file contains a list of pre-defined comparisons between groups (e.g. treatment versus controls, mutant versus wild-type). If you do not have a comparisons file, you can create comparisons using the interactive comparison builder.", closable = FALSE),
                style = "align-items: end"
              ),
              bslib::card(
                fileInputArea(
                  ns("contrasts_csv"),
                  div(
                    style = "display: flex; flex-direction: column; gap: 10px;",
                    shiny::h4("Upload comparisons.csv (optional)", class = "mb-0")
                  ),
                  multiple = FALSE,
                  accept = c(".csv"),
                  width = "100%"
                ),
                style = "background-color: aliceblue; border: 0.07rem dashed steelblue;"
              ),
              action_buttons2
            )
          } else {
            bslib::layout_columns(
              col_widths = c(8, 4),
              TableModuleUI(
                ns("contrasts_datasets"),
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
                  plotOutput(ns("contrasts_stats")),
                  br(), hr(),
                  "Summary:",
                  br(),
                  check_to_html(
                    checklist$contrasts.csv$checks,
                    pass_msg = "All contrasts checks passed",
                    null_msg = "Contrasts checks not run yet.
                                  Fix any errors with contrasts first.",
                    details = FALSE
                  ),
                  check_to_html(checklist$samples_contrasts$checks,
                    pass_msg = "All contrasts-samples checks passed",
                    null_msg = "Contrasts-samples checks not run yet.
                              Fix any errors with contrasts or samples first.",
                    details = FALSE
                  )
                  ## preview_module_legend
                )
              ),
              action_buttons2
            )
          }
        )

      div(
        # if run_build_comparisons is clicked, then show the contrasts
        bslib::as_fill_carrier(),
        if (show_comparison_builder()) {
          bslib::layout_columns(
            col_widths = 12,
            ## height = "calc(100vh - 340px)",
            heights_equal = "row",
            bs_alert(HTML("To <b>create comparisons</b>, choose a phenotype, then create groups by dragging conditions to the 'Main group' or 'Control group' boxes, give a name and click 'add to list'. You can also try 'auto-detect comparisons'. If you have a file with pre-defined comparisons, you can upload this below.")),
            upload_module_makecontrast_ui(ns("makecontrast")),
            action_buttons1
          )
        },
        if (!show_comparison_builder()) {
          comparison_builder_ui
        }
      )
      # }
    })

    output$contrasts_stats <- renderPlot({
      ct <- uploaded$contrasts.csv
      shiny::req(nrow(ct))
      plotPhenoDistribution(data.frame(ct))
    })

    # error pop-up alert
    observeEvent(
      {
        list(
          checklist$contrasts.csv$checks,
          checklist$samples_contrasts$checks
        )
      },
      {
        checks1 <- checklist$samples.csv$checks
        checks2 <- checklist$samples_contrasts$checks

        if (length(checks1) > 0 || length(checks2) > 0) {
          err.html <- ""
          if (length(checks1) > 0) {
            err1 <- check_to_html(
              checks1,
              pass_msg = tspan("All contrasts checks passed", js = FALSE),
              null_msg = tspan("Contrasts checks not run yet.
                            Fix any errors with counts first.", js = FALSE),
              details = TRUE
            )
            err.html <- paste(err.html, err1)
          }
          if (length(checks2) > 0) {
            err2 <- check_to_html(
              checks2,
              pass_msg = tspan("All samples-contrasts checks passed", js = FALSE),
              null_msg = tspan("Samples-contrasts checks not run yet.
                        Fix any errors with samples or contrasts first.", js = FALSE),
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


    # control state of comparison builder
    observeEvent(input$run_build_comparisons, {
      show_comparison_builder(TRUE)
    })

    observeEvent(input$goUploadComparison, {
      selected_contrast_input(TRUE)
      show_comparison_builder(FALSE)
    })

    observeEvent(input$goOnlineComparison, {
      selected_contrast_input(TRUE)
      show_comparison_builder(TRUE)
    })

    # pass counts to uploaded when uploaded
    observeEvent(input$contrasts_csv, {
      # check if contrasts is csv (necessary due to drag and drop of any file)
      ext <- tools::file_ext(input$contrasts_csv$name)[1]
      if (ext != "csv") {
        shinyalert::shinyalert(
          title = "File format not supported.",
          text = "Please make sure the file is a CSV file.",
          type = "error"
        )
        return()
      }

      # if contrasts not in file name, give warning and return
      if (!grepl("contrast", input$contrasts_csv$name, ignore.case = TRUE) && !grepl("comparison", input$contrasts_csv$name, ignore.case = TRUE)) {
        shinyalert::shinyalert(
          title = "Comparison not in filename.",
          text = "Please make sure the file name contains 'comparison', such as comparison_dataset.csv or comparison.csv.",
          type = "error"
        )
        return()
      }

      ct <- tryCatch(
        {
          playbase::read.as_matrix(input$contrasts_csv$datapath)
        },
        error = function(w) {
          NULL
        }
      )
      if (is.null(ct)) {
        data_error_modal(
          path = input$contrasts_csv$datapath,
          data_type = "contrasts"
        )
      } else {
        ## IK: should we do contrasts.convertToLabelMatrix here
        ## already?  to allow for short/old formats?
        ## samples <- checked_samples()$SAMPLES
        samples <- uploaded$samples.csv
        new.ct <- try(playbase::contrasts.convertToLabelMatrix(
          contrasts = ct, samples = samples
        ))
        if (!"try-error" %in% class(new.ct)) {
          ct <- new.ct
        }
        uploaded$contrasts.csv <- ct
      }
    })

    observeEvent(input$remove_contrasts, {
      uploaded$contrasts.csv <- NULL
      checklist$samples_contrasts$checks <- NULL
    })

    observeEvent(input$load_example, {
      uploaded$contrasts.csv <- playbase::CONTRASTS
    })

    modified_ct <- upload_module_makecontrast_server(
      id = "makecontrast",
      phenoRT = reactive(checked_samples()$SAMPLES),
      contrRT = reactive(checked_contrasts()$matrix),
      countsRT = reactive(checked_counts()$COUNTS),
      upload_wizard = upload_wizard,
      show_comparison_builder = show_comparison_builder,
      autocontrast = reactive(input$autocontrast)
    )

    TableModuleServer(
      "contrasts_datasets",
      func = table.RENDER,
      selector = "none"
    )

    return(modified_ct)
  })
}

# convert list of checks to html tags for display in the data preview modal
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
