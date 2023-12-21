# NOTES
# SAMPLES
#- all tabs should be the same size
#- convert column types from numeric to categorical
#- remove column
#- convert columns to categorical via binning and cuts
#
# CONTRASTS
#- make contrasts table editable
#
#- Save updated files
#
#- put the errors in the same tab
#- make the modal the whole screen and put the errors on the side
#

# no ui right now because the preview is in a modal dialog
upload_module_preview_ui <- function(id) {
  ns <- NS(id)
  tagList()
}

upload_module_preview_server <- function(id, uploaded, checklist, checkTables) {
  moduleServer(
    id,
    function(input, output, session) {
      
      # ever/y time something is uploaded, it can be previewed
      observeEvent( c(uploaded$last_uploaded), {
        checkTables()
        has_counts <- !is.null(uploaded$counts.csv)
        has_samples <- !is.null(uploaded$samples.csv)
        has_contrasts <- !is.null(uploaded$contrasts.csv)
        
        tabs <- list(id = session$ns("preview_panel"))
        if (has_counts) {
          tabs <- c(
            tabs,
            list(tabPanel(
              "Counts",
              fluidRow(
                column(
                  width = 9,
                  upload_table_preview_counts_ui(session$ns("counts_preview"),
                    height = c("100%", TABLE_HEIGHT_MODAL),
                    width = c("auto", "100%"),
                    title = "Uploaded Counts",
                    info.text = "This is the uploaded counts data.",
                    caption = "This is the uploaded counts data."
                    )
                ),
                column(
                  width = 3,
                  class = "px-4",
                  "Summary:", br(),
                  check_to_html(checklist$counts.csv$checks,
                    pass_msg = "All counts checks passed",
                    null_msg = "Counts checks not run yet.
                                Fix any errors with counts first."
                  ),
                  ifelse(!has_samples, "",
                    check_to_html(checklist$samples_counts$checks,
                        pass_msg = "All samples-counts checks passed",
                        null_msg = "Samples-counts checks not run yet.
                                Fix any errors with samples or counts first."
                      )
                  )
                )
              )
            ))
          )
        }
        if (has_samples) {
          tabs <- c(
            tabs,
            list(tabPanel(
              "Samples",
              fluidRow(
                column(
                  width = 9,
                  upload_table_preview_samples_ui(session$ns("samples_preview"),
                    height = c("100%", TABLE_HEIGHT_MODAL),
                    width = c("auto", "100%"),
                    title = "Uploaded Samples",
                    info.text = "This is the uploaded samples data.",
                    caption = "This is the uploaded samples data."
                  )
                ),
                column(
                  width = 3,
                  class = "px-4",
                  "Summary:", br(),
                  check_to_html(checklist$samples.csv$checks,
                    pass_msg = "All samples checks passed",
                    null_msg = "Samples checks not run. Fix any
                                errors with samples first."
                  ),
                  ifelse(!has_counts, "",
                    check_to_html(
                      check = checklist$samples_counts$checks,
                      pass_msg = "All samples-counts checks passed",
                      null_msg = "Samples-counts checks not run yet.
                                  Fix any errors with samples or counts first."
                    )
                  ),
                  ifelse(!has_contrasts, "",
                    check_to_html(checklist$samples_contrasts$checks,
                      pass_msg = "All samples-contrasts checks passed",
                      null_msg = "Samples-contrasts checks not run yet.
                                Fix any errors with samples or contrasts first."
                    )
                  )
                )
              )
            ))
          )
        }
        if (has_contrasts) {
          tabs <- c(
            tabs,
            list(tabPanel(
              "Contrasts",
              fluidRow(
                column(
                  width = 9,
                  upload_table_preview_contrasts_ui(session$ns("contrasts_preview"),
                    height = c("100%", TABLE_HEIGHT_MODAL),
                    width = c("auto", "100%"),
                    title = "Uploaded Contrasts",
                    info.text = "This is the uploaded contrasts data.",
                    caption = "This is the uploaded contrasts data."
                  )
                ),
                column(
                  width = 3,
                  class = "px-4",
                  "Summary:", br(),
                  check_to_html(checklist$contrasts.csv$checks,
                    pass_msg = "All contrasts checks passed",
                    null_msg = "Contrasts checks not run. Fix any errors
                                with contrasts first."
                  ),
                  ifelse(!has_samples, "",
                    check_to_html(checklist$samples_contrasts$checks,
                      pass_msg = "All samples-contrasts checks passed",
                      null_msg = "Samples-contrasts checks not run yet.
                                  Fix any errors with samples or contrasts first."
                    )
                  )
                )
              )
            ))
          )
        }

        shiny::showModal(
          shiny::modalDialog(
            title = "Data Upload Preview",
            label = "Data Upload Preview",
            do.call(tabsetPanel, tabs),
            fluidRow(
              column(
                12,
                class = "px-4",
                span(style = "color: orange", "Orange"),
                span("= warning but data will still be uploaded. "),
                span(style = "color:red", "Red"),
                span("= error and data will not be uploaded.")
              )
            ),
            footer = div(
              style = "float: right",
              shiny::actionButton(
                session$ns("cancel_upload"),
                label = "Cancel",
                style = "display: inline-block; margin-left: 15px;",
                class = "btn-danger"
              ),
              shiny::actionButton(
                session$ns("ok_upload"),
                label = "OK",
                style = "display: inline-block",
                class = "btn-info"
              )
            ),
            easyClose = FALSE,
            size = "xl"
          ) %>%
            tagAppendAttributes(
              style = "min-height: 90%; min-width: 90%",
              .cssSelector = ".modal-dialog"
            )
        )
      },
      ignoreNULL = TRUE
      )

      upload_table_preview_counts_server("counts_preview",
        uploaded,
        scrollY = "calc(35vh - 140px)"
      )
      upload_table_preview_samples_server("samples_preview",
        uploaded,
        scrollY = "calc(35vh - 140px)"
      )
      upload_table_preview_contrasts_server("contrasts_preview",
        uploaded,
        scrollY = "calc(35vh - 140px)"
      )

      observeEvent(input$cancel_upload,
      {
        shiny::removeModal()
        has_counts <- "counts.csv" %in% uploaded$last_uploaded
        has_samples <- "samples.csv" %in% uploaded$last_uploaded
        has_contrasts <- ("contrasts.csv" %in% uploaded$last_uploaded) &
          (!is.null(uploaded$contrasts.csv))
        uploaded[["counts.csv"]] <- NULL
        uploaded[["samples.csv"]] <- NULL
        uploaded[["contrasts.csv"]] <- NULL
        uploaded[["pgx"]] <- NULL
        uploaded[["last_uploaded"]] <- NULL
        uploaded[["checklist"]] <- NULL
        checklist[["counts.csv"]] <- NULL
        checklist[["samples.csv"]] <- NULL
        checklist[["contrasts.csv"]] <- NULL
        checklist[["samples_counts"]] <- NULL
        checklist[["samples_contrasts"]] <- NULL
      },
      ignoreInit = TRUE
      )

      observeEvent(input$ok_upload,
      {
        shiny::removeModal()
      },
      ignoreInit = TRUE
      )
    }
  )
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


upload_table_preview_counts_ui <- function(
    id,
    width,
    height,
    title,
    info.text,
    caption) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    width = width,
    height = height,
    title = title,
    info.text = info.text,
    caption = caption,
    label = "",
    show.maximize = FALSE
  )
}

upload_table_preview_counts_server <- function(id,
                                               uploaded,
                                               scrollY) {
  moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      shiny::req(uploaded$counts.csv)
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
      req(dt)

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
        DT::formatRound( columns = 1:ncol(dt), digits=3) %>%     
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    TableModuleServer(
      "datasets",
      func = table.RENDER,
      selector = "none"
    )
  }) ## end of moduleServer
} ## end of server

upload_table_preview_samples_ui <- function(
    id,
    width,
    height,
    title,
    info.text,
    caption) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    width = width,
    height = height,
    title = title,
    info.text = info.text,
    caption = caption,
    label = "",
    show.maximize = FALSE
  )
}

upload_table_preview_samples_server <- function(id,
                                                uploaded,
                                                scrollY) {
  moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      shiny::req(uploaded$samples.csv)
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
      dt
    })

    table.RENDER <- function() {
      dt <- table_data()
      req(dt)
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

    TableModuleServer(
      "datasets",
      func = table.RENDER,
      selector = "none"
    )
  }) ## end of moduleServer
} ## end of server

upload_table_preview_contrasts_ui <- function(
    id,
    width,
    height,
    title,
    info.text,
    caption) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    width = width,
    height = height,
    title = title,
    info.text = info.text,
    caption = caption,
    label = "",
    show.maximize = FALSE
  )
}

upload_table_preview_contrasts_server <- function(id,
                                                  uploaded,
                                                  scrollY) {
  moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      shiny::req(uploaded$contrasts.csv)
      dt <- uploaded$contrasts.csv |> data.frame(check.names = FALSE)
      if (NCOL(dt) == 0) dt <- cbind(dt, " " = NA)
      dt
    })

    table.RENDER <- function() {
      dt <- table_data()
      req(dt)
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

    TableModuleServer(
      "datasets",
      func = table.RENDER,
      selector = "none"
    )
  }) ## end of moduleServer
} ## end of server


# function to format datatable for highlighting and selecting columns instead of rows
# not currently used
make_preview_table <- function(tbl) {
  # callback for highlighting column instead of row
  js <- c(
    "table.on('mouseenter', 'td', function () {
                        // Remove highlight from all columns
                        table
                        .columns()
                        .nodes()
                        .flatten()  // Reduce to a 1D array
                        .to$()      // Convert to a jQuery object
                        .removeClass( 'highlight' );

                                // Add highlight to mouseover column
                        table
                        .column( this )
                        .nodes()
                        .to$()      // Convert to a jQuery object
                        .addClass( 'highlight' );
                    });",
    "table.on('mouseleave', 'td', function () {
                        // Remove highlight from all columns
                        table
                        .columns()
                        .nodes()
                        .flatten()  // Reduce to a 1D array
                        .to$()      // Convert to a jQuery object
                        .removeClass( 'highlight' );
                    });"
  )

  DT::datatable(tbl,
    class = "compact",
    rownames = TRUE,
    options = list(
      dom = "rtp",
      pageLength = 20
    ),
    callback = DT::JS(js),
    selection = list(target = "column", mode = "single")
  ) %>%
    DT::formatStyle(0,
      target = "row",
      fontSize = "12px", lineHeight = "70%"
    )
}
