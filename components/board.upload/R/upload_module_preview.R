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
  caption,
  checked_counts,
  checked_samples,
  cjecked
  ) 
  {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

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

    output$table_counts <- shiny::renderUI(
      div(
        bslib::as_fill_carrier(),
        div(
          style = "display: flex; justify-content: space-between; margin-bottom: 20px;",
          div(
            if(!is.null(uploaded$counts.csv)){
              shiny::actionButton(
                ns("remove_counts"),
                "Remove input",
                icon = icon("trash-can"),
                class = "btn btn-danger"
              )
            }
        ),
        div(
          actionButton(
            ns("load_example"), "Load Example",
            class = "btn btn-info"
            ),
          shiny::actionButton(
            ns("check_documentation_counts"),
            "Check Documentation",
            class = "btn btn-primary",
            onclick ="window.open('https://omicsplayground.readthedocs.io/en/latest/dataprep/counts/', '_blank')"
            )
        )
        ),
        if(is.null(uploaded$counts.csv)){
        bslib::layout_columns(
          bslib::card(
            fileInputArea(
              ns("counts_csv"),
              shiny::h4("Choose counts.csv", class='mb-0'),
              multiple = FALSE,
              accept = c(".csv"),
              width = "800px" #TODO this is a hack to make the file input area wider
            )
          )
        )
      }else{
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
            legend
            )
          )
      )
      }
      )
    )

    # pass counts to uploaded when uploaded
    observeEvent(input$counts_csv, {

      # if counts not in file name, give warning and return
      if(!grepl("count", input$counts_csv$name, ignore.case = TRUE) && !grepl("expression", input$counts_csv$name, ignore.case = TRUE)){
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

      #TODO this needs review!
      # if(IS_EXPRESSION){
      #   df <- 2**df0 
      #   if(min(df0,na.rm=TRUE) > 0) df <- df - 1
      # }

      # #TODO this needs review!
      # if (IS_COUNT){
      #   check.log <- ( all( df0 < 60) || min(df0, na.rm=TRUE) < 0 )
      #   if(check.log) {
      #     shinyalert::shinyalert(
      #       title = "",
      #       text = "Count matrix provided but log2-values detected. Converting log2 value to intensities.",
      #       type = "warning"
      #       )
      #     #TODO this conversion should be optional ,upon confirmation of user
      #     df <- 2**df0 - 1
      #     if(min(df0, na.rm=TRUE) > 0) df <- df - 1
      #   }
      # }
      uploaded$counts.csv <- df

    })

    observeEvent(input$remove_counts, {

      delete_all_files_counts <- function(value) {
        if(input$alert_delete_counts){
          uploaded$counts.csv <- NULL
          uploaded$samples.csv <- NULL
          uploaded$contrasts.csv <- NULL
        }
      }
      
      # if samples is not null, warn user that it will be deleted
      if(!is.null(uploaded$samples.csv) || !is.null(uploaded$contrasts.csv)){
        shinyalert::shinyalert(
          inputId = "alert_delete_counts",
          title = "Warning",
          text = "Removing counts will also remove samples and contrasts. Do you want to proceed?",
          type = "warning",
          showCancelButton = TRUE,
          closeOnEsc = FALSE,
          callbackR = delete_all_files_counts,
          confirmButtonText = "Yes, remove all files.",
          cancelButtonText = "No, cancel deletion."
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

upload_table_preview_samples_ui <- function(id) {
  ns <- shiny::NS(id)
  uiOutput(ns("table_samples"), fill=TRUE)
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
  caption
  ) 
  {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

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
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    output$table_samples <- shiny::renderUI(
      div(
        bslib::as_fill_carrier(),
        div(
          style = "display: flex; justify-content: space-between; margin-bottom: 20px;",
          div(
            if(!is.null(uploaded$samples.csv)){
              shiny::actionButton(
                ns("remove_samples"),
                "Remove input",
                icon = icon("trash-can"),
                class = "btn btn-danger"
              )
            }
        ),
        div(
          actionButton(
            ns("load_example"), "Load Example",
            class = "btn btn-info"
            ),
          actionButton(
            ns("check_documentation_samples"),
            "Check Documentation",
            class = "btn btn-primary",
            onclick ="window.open('https://omicsplayground.readthedocs.io/en/latest/dataprep/samples/', '_blank')"
            )
        )
        ),
        if(is.null(uploaded$samples.csv)){
        bslib::layout_columns(
          bslib::card(
            fileInputArea(
              ns("samples_csv"),
              shiny::h4("Choose samples.csv", class='mb-0'),
              multiple = FALSE,
              accept = c(".csv"),
              width = "800px" #TODO this is a hack to make the file input area wider
            )
          )
        )
      }else{
         bslib::layout_columns(
        col_widths = c(9, 3),
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
          div(
            "Summary:",
            br(),
            check_to_html(
              checklist$samples.csv$checks,
                pass_msg = "All samples checks passed",
                null_msg = "Samples checks not run yet.
                            Fix any errors with samples first."
              ),
            check_to_html(checklist$samples_counts$checks,
                pass_msg = "All samples-counts checks passed",
                null_msg = "Samples-counts checks not run yet.
                        Fix any errors with samples or counts first."
              ),
            legend
            )
          )
      )
      }
      )
    )

    # pass counts to uploaded when uploaded
    observeEvent(input$samples_csv, {
      # if samples not in file name, give warning and return
      if(!grepl("sample", input$samples_csv$name, ignore.case = TRUE)){
        shinyalert::shinyalert(
          title = "Samples not in filename.",
          text = "Please make sure the file name contains 'samples', such as samples_dataset.csv or samples.csv.",
          type = "error"
        )
        return()
      }
      
      uploaded$samples.csv <- playbase::read.as_matrix(input$samples_csv$datapath)

    })

    observeEvent(input$remove_samples, {
      
      delete_all_files_samples <- function(value) {
        if(input$alert_delete_samples){
          uploaded$samples.csv <- NULL
          uploaded$contrasts.csv <- NULL
        }
      }

      # if samples is not null, warn user that it will be deleted
      if(!is.null(uploaded$contrasts.csv)){
        
        shinyalert::shinyalert(
          inputId = "alert_delete_samples",
          title = "Warning",
          text = "Removing samples will also remove contrasts. Do you want to proceed?",
          type = "warning",
          showCancelButton = TRUE,
          closeOnEsc = FALSE,
          callbackR = delete_all_files_samples,
          confirmButtonText = "Yes, remove samples and contrasts.",
          cancelButtonText = "No, cancel deletion."
        )

        } else {
          uploaded$samples.csv <- NULL
        }
    })

    observeEvent(input$load_example, {
      uploaded$samples.csv <- playbase::SAMPLES
    })

    TableModuleServer(
      "samples_datasets",
      func = table.RENDER,
      selector = "none"
    )
  })
}

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
  upload_wizard
  )
  {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    table_data <- shiny::reactive({
      shiny::req(uploaded$contrasts.csv)
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

    output$table_contrasts <- shiny::renderUI(
      div(
        # if run_build_comparisons is clicked, then show the contrasts
        bslib::as_fill_carrier(),
        if(show_comparison_builder()){
          bslib::layout_columns(
          col_widths = 12,
          # height = "calc(100vh - 340px)",
          heights_equal = "row",
          upload_module_makecontrast_ui(ns("makecontrast")),
          #bs_alert(HTML("Here, you can interactively <b>create comparisons</b> (also called 'contrasts'). Choose a phenotype, then create groups by dragging conditions to the boxes of the 'main' or 'control' group. Give the contrast a name (please keep it short!) and then click 'add comparison'. If you are feeling lucky, you can also try 'auto-comparisons'."))
          )
        } else {
          div(
          div(
            style = "display: flex; justify-content: space-between; margin-bottom: 20px;",
            div(
            if(!is.null(uploaded$contrasts.csv)){
                shiny::actionButton(
                  ns("remove_contrasts"),
                  "Remove input",
                  icon = icon("trash-can"),
                  class = "btn btn-danger"
                )
              },
              actionButton(
                  ns("run_build_comparisons"), "Run Comparison builder",
                  class = "btn btn-warning"
                ),
              ),
            div(
            actionButton(
              ns("load_example"), "Load Example",
              class = "btn btn-info"
              ),
            actionButton(
              ns("check_documentation_contrasts"),
              "Check Documentation",
              class = "btn btn-primary",
              onclick ="window.open('https://omicsplayground.readthedocs.io/en/latest/dataprep/contrasts/', '_blank')"
              )
          )),
          if(is.null(uploaded$contrasts.csv)){
          bslib::layout_columns(
            bslib::card(
              fileInputArea(
                ns("contrasts_csv"),
                shiny::h4("Choose contrasts.csv", class='mb-0'),
                multiple = FALSE,
                accept = c(".csv"),
                width = "800px" #TODO this is a hack to make the file input area wider
              )
            )
          )
        }else{
          bslib::layout_columns(
          col_widths = c(9, 3),
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
              "Summary:",
              br(),
              check_to_html(
                checklist$contrasts.csv$checks,
                  pass_msg = "All contrasts checks passed",
                  null_msg = "Contrasts checks not run yet.
                              Fix any errors with contrasts first."
                ),
              check_to_html(checklist$samples_contrasts$checks,
                  pass_msg = "All contrasts-samples checks passed",
                  null_msg = "Contrasts-samples checks not run yet.
                          Fix any errors with contrasts or samples first."
                ),
              legend
              )
            )
        )
        }
        )
        }
    ))

    # control state of comparison builder
    observeEvent(input$run_build_comparisons, {
      show_comparison_builder(TRUE)
    })

    # pass counts to uploaded when uploaded
    observeEvent(input$contrasts_csv, {
      
      # if contrasts not in file name, give warning and return
      if(!grepl("contrast", input$contrasts_csv$name, ignore.case = TRUE)){
        shinyalert::shinyalert(
          title = "Contrasts not in filename.",
          text = "Please make sure the file name contains 'contrasts', such as contrasts_dataset.csv or contrasts.csv.",
          type = "error"
        )
        return()
      }
      
      uploaded$contrasts.csv <- playbase::read.as_matrix(input$contrasts_csv$datapath)

    })

    observeEvent(input$remove_contrasts, {
      uploaded$contrasts.csv <- NULL
    })

    observeEvent(input$load_example, {
      uploaded$contrasts.csv <- playbase::CONTRASTS
    })

    modified_ct <- upload_module_makecontrast_server(
      id = "makecontrast",
      phenoRT = reactive(checked_samples()$matrix),
      contrRT = reactive(checked_contrasts()$matrix),
      countsRT = reactive(checked_counts()$matrix),
      upload_wizard = upload_wizard,
      show_comparison_builder = show_comparison_builder
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


legend <-  shiny::div(
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


fileInputArea <- function(inputId, label, multiple = FALSE, accept = NULL,
                          width = NULL, buttonLabel = "Browse...", placeholder = "No file selected") {
  restoredValue <- restoreInput(id = inputId, default = NULL)

  # Catch potential edge case - ensure that it's either NULL or a data frame.
  if (!is.null(restoredValue) && !is.data.frame(restoredValue)) {
    warning("Restored value for ", inputId, " has incorrect format.")
    restoredValue <- NULL
  }

  if (!is.null(restoredValue)) {
    restoredValue <- toJSON(restoredValue, strict_atomic = FALSE)
  }

  inputTag <- tags$input(
    id = inputId,
    name = inputId,
    type = "file",
    # Don't use "display: none;" style, which causes keyboard accessibility issue; instead use the following workaround: https://css-tricks.com/places-its-tempting-to-use-display-none-but-dont/
    style = "position: absolute !important; top: -99999px !important; left: -99999px !important;",
    `data-restore` = restoredValue
  )

  if (multiple) {
    inputTag$attribs$multiple <- "multiple"
  }
  if (length(accept) > 0) {
    inputTag$attribs$accept <- paste(accept, collapse = ",")
  }

  div(
    class = "form-group", #shiny-input-container w-100",
    style = htmltools::css(width = htmltools::validateCssUnit(width)),
    shiny:::shinyInputLabel(inputId, ""),
    div(
      class = "input-group mb-3",
      # input-group-prepend is for bootstrap 4 compat
      tags$label(
        class = "input-group-btn input-group-prepend w-100",
        span(
          class = "btn btn-area w-100", inputTag,
          div(tags$image(src = icon_encoded, width = "80px;"), style = "margin-top: 2rem;"),
          div(p(label), style = "font-size: 1.2rem; font-weight: 700; padding-top: 2rem;"),
          div(p(buttonLabel), style = "font-size: 1rem; font-weight: 400; margin-bottom: 2rem;")
        )
      )
    ),
    tags$div(
      id = paste(inputId, "_progress", sep = ""),
      class = "progress active shiny-file-input-progress",
      tags$div(class = "progress-bar")
    )
  )
}

# Icon from <https://icons.getbootstrap.com/icons/upload/>
icon_file <- tempfile(fileext = ".svg")
writeLines('
<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="#495057" class="bi bi-upload" viewBox="0 0 16 16">
  <path d="M.5 9.9a.5.5 0 0 1 .5.5v2.5a1 1 0 0 0 1 1h12a1 1 0 0 0 1-1v-2.5a.5.5 0 0 1 1 0v2.5a2 2 0 0 1-2 2H2a2 2 0 0 1-2-2v-2.5a.5.5 0 0 1 .5-.5z"/>
  <path d="M7.646 1.146a.5.5 0 0 1 .708 0l3 3a.5.5 0 0 1-.708.708L8.5 2.707V11.5a.5.5 0 0 1-1 0V2.707L5.354 4.854a.5.5 0 1 1-.708-.708l3-3z"/>
</svg>',
  con = icon_file
)

icon_encoded <- xfun::base64_uri(icon_file)