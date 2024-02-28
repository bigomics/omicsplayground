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
      if(is.null(input$counts_csv)){
        bslib::layout_columns(
          bslib::card(
            fileInput2(
              ns("counts_csv"),
              shiny::h4("Choose counts.csv", class='mb-0'),
              multiple = FALSE,
              buttonClass = "btn-primary",
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
            check_to_html(
              checklist$samples_counts$checks,
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