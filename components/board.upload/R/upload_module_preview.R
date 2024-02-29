upload_table_preview_counts_ui <- function(id) {
  
  ns <- shiny::NS(id)
  uiOutput(ns("table_counts"))
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
  caption
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
        div(
          style = "display: flex; justify-content: space-between;",
          div(
            shiny::actionButton(
              ns("remove_counts"),
              "Remove input", 
              icon = icon("trash-can"),
              class = "btn btn-outline-danger"
              )
        ),
        div(
          actionButton(
            ns("load_example"), "Load Example",
            class = "btn btn-outline-info"
            ),
          actionButton(
            ns("check_documentation"),
            "Check Documentation",
            class = "btn btn-outline-primary"
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
              accept = c(".csv")
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
      
      
      df0 <- playbase::read.as_matrix(input$counts_csv$datapath)
      df <- df0
      # file.copy(df0, file.path(raw_dir(), "raw_counts.csv"))

      IS_EXPRESSION <- grepl("expression", input$counts_csv$datapath, ignore.case = TRUE)
      IS_COUNT <- grepl("count", input$counts_csv$datapath, ignore.case = TRUE)

      
      
      #TODO this needs review!
      if(IS_EXPRESSION){
        df <- 2**df0 
        if(min(df0,na.rm=TRUE) > 0) df <- df - 1
      }

      #TODO this needs review!
      if (IS_COUNT){
        check.log <- ( all( df0 < 60) || min(df0, na.rm=TRUE) < 0 )
        if(check.log) {
          shinyalert::shinyalert(
            title = "",
            text = "Count matrix provided but log2-values detected. Converting log2 value to intensities.",
            type = "warning"
            )
          df <- 2**df0 - 1
          if(min(df0, na.rm=TRUE) > 0) df <- df - 1
        }
      }
      uploaded$counts.csv <- df

    })

    observeEvent(input$remove_counts, {
      uploaded$counts.csv <- NULL
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

upload_table_preview_samples_ui <- function(
    id,
    width,
    height,
    title,
    info.text,
    caption) {
  ns <- shiny::NS(id)

  bslib::layout_columns(
    col_widths = c(9, 3),
    TableModuleUI(
      ns("datasets"),
      width = width,
      height = height,
      title = title,
      info.text = info.text,
      caption = caption,
      label = "",
      show.maximize = FALSE
    ),
    bslib::card(
      uiOutput(ns("checklist"))
    )
  )
}

upload_table_preview_samples_server <- function(id,
                                                uploaded,
                                                checklist,
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

    output$checklist <- renderUI({
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
    })

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


  bslib::layout_columns(
    col_widths = c(9, 3),
    TableModuleUI(
      ns("datasets"),
      width = width,
      height = height,
      title = title,
      info.text = info.text,
      caption = caption,
      label = "",
      show.maximize = FALSE
    ),
    bslib::card(
      uiOutput(ns("checklist"))
    )
  )
}

upload_table_preview_contrasts_server <- function(id,
                                                  uploaded,
                                                  checklist,
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


    output$checklist <- renderUI({
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
              pass_msg = "All samples-contrasts checks passed",
              null_msg = "Samples-contrasts checks not run yet.
                      Fix any errors with samples or contrasts first."
          ),
          legend
        )
    })

    TableModuleServer(
      "datasets",
      func = table.RENDER,
      selector = "none"
    )
  }) ## end of moduleServer
} ## end of server


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