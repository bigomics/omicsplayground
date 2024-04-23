##---------------------------------------------------
## COUNTS
##---------------------------------------------------

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

    output$table_counts <- shiny::renderUI({

      action_buttons <- div(
        style = "display: flex; justify-content: left; margin: 8px;",
        if(is.null(uploaded$counts.csv)){
          div(               
            actionButton(
              ns("load_example"), "Load example data",
              class = "btn-sm btn-outline-primary m-1"
            )
          )
        } else {
          shiny::actionButton(
            ns("remove_counts"),
            "Cancel",
            icon = icon("trash-can"),
            class = "btn-sm btn-outline-danger m-1"
            )
        },
        shiny::actionButton(
          ns("check_documentation_counts"),
          "Read documentation",
          class = "btn-sm btn-outline-primary m-1",
          onclick ="window.open('https://omicsplayground.readthedocs.io/en/latest/dataprep/counts/', '_blank')"
          )
      )

      alert <- bs_alert("The counts file (counts.csv) contains the measurements (genes, proteins, etc..) for all samples. The file should be a tabular text file (.csv), where each row corresponds to a feature (i.e. genes or proteins) and each column corresponds to a sample.", closable = FALSE) 
      htmltools::tagAppendChild(div(), alert)
      
      div(
        bslib::as_fill_carrier(),
        ## style = "width: 100%; display: flex; justify-content: space-between; margin-bottom: 8px;",
        style = "width: 100%; display: flex; ",
        if(is.null(uploaded$counts.csv)){
          bslib::layout_columns(          
            col_widths = c(-3,6,-3),
            row_heights = list("auto","auto","auto"),
            gap = "0.5rem",
            alert,
            bslib::card( 
              fileInputArea(
                ns("counts_csv"),
                shiny::h4("Upload counts.csv", class='mb-0'),
                multiple = FALSE,
                accept = c(".csv"),
                width = "100%"
              ),
              style = "background-color: aliceblue; border: 0.1rem dashed steelblue;"
            ),
            action_buttons            
          )
        } else {
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
            ),
            action_buttons          
          )
        }, ## end of if-else
      ) ## end of div

    })

    observeEvent(input$upload_counts, {
      shinyjs::click("counts_csv")
    })
    
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

##---------------------------------------------------
## SAMPLES
##---------------------------------------------------

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

    output$table_samples <- shiny::renderUI({

      action_buttons <- div(
        style = "display: flex; justify-content: left; margin-bottom: 8px;",
        div(
          if(!is.null(uploaded$samples.csv)){
            shiny::actionButton(
              ns("remove_samples"),
              "Cancel",
              icon = icon("trash-can"),
              class = "btn-sm btn-outline-danger m-1"
              )
          } else {
            shiny::actionButton(
              ns("load_example"), "Load example data",
              class = "btn-sm btn-outline-primary m-1"
            )
          }
        ),
        div(
          shiny::actionButton(
            ns("check_documentation_samples"),
            "Read documentation",
              class = "btn-sm btn-outline-primary m-1",
            onclick ="window.open('https://omicsplayground.readthedocs.io/en/latest/dataprep/samples/', '_blank')"
          )
        )
      )
      
      div(
        bslib::as_fill_carrier(),
        if(is.null(uploaded$samples.csv)) {
          bslib::layout_columns(
            col_widths = c(-3, 6, -3),                   
            gap = "0.3rem",
            div( bs_alert("The samples file (samples.csv) contains the phenotypic information of your samples. The file should be a tabular text file (csv) with the samples in the rows and the phenotypic data (metadata) in the columns. The first column contains the sample names, which must be unique, and has to match the names given in the header of the read counts file.", closable = FALSE), style = "margin-bottom: -50px;"),
            bslib::card(
              fileInputArea(
                ns("samples_csv"),
                shiny::h4("Choose samples.csv", class='mb-0'),
                multiple = FALSE,
                accept = c(".csv"),
                width = "100%"
              ),
              style = "background-color: aliceblue; border: 0.1rem dashed steelblue;"
            ),
            action_buttons
          )
        } else {
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
            ),
            action_buttons            
          )
        }
      )
    })

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
          confirmButtonText = "Yes, remove samples and contrasts",
          cancelButtonText = "Cancel"
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

##---------------------------------------------------
## CONTRASTS
##---------------------------------------------------

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

    output$table_contrasts <- shiny::renderUI({
      
      if(selected_contrast_input() == FALSE) {
        # ask user if preferrence is upload contrast or create contrast online
        div(
          style = "display: flex; gap: 80px; flex-direction: column; justify-content: center; align-items: center; height: 50%;",
          div(h4("Please choose one of the following options:")),
          div(
            style = "display: flex; justify-content: center; gap: 20px;",
            actionButton(
            ns("goUploadComparison"),
            label = " Upload comparisons file",
            class = "btn btn-primary m-1",
            icon("upload")
            ),
          actionButton(
            ns("goOnlineComparison"),
            label = " Build my comparisons",
            class = "btn btn-primary m-1",
            icon("pen-to-square")
            )
          )
          
        )
      } else {
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
            onclick ="window.open('https://omicsplayground.readthedocs.io/en/latest/dataprep/contrasts/', '_blank')"
          )
        )

        action_buttons2 <- div(
          style = "display: flex; justify-content: left; margin-bottom: 8px;",
          if(!is.null(uploaded$contrasts.csv)){
            shiny::actionButton(
              ns("remove_contrasts"),
              "Cancel",
              icon = icon("trash-can"),
              class = "btn-sm btn-outline-danger m-1"
              )
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
            onclick ="window.open('https://omicsplayground.readthedocs.io/en/latest/dataprep/contrasts/', '_blank')"
              )
        )
        
        
        div(
          # if run_build_comparisons is clicked, then show the contrasts
          bslib::as_fill_carrier(),
          if(show_comparison_builder()){
            bslib::layout_columns(
              col_widths = 12,
              ## height = "calc(100vh - 340px)",
              heights_equal = "row",
              bs_alert(HTML("To <b>create comparisons</b>, choose a phenotype, then create groups by dragging conditions to the 'Main' or 'Control' group. Give the contrast a name and click 'add'. You can also try 'auto-detect comparisons'.")),
              upload_module_makecontrast_ui(ns("makecontrast")),
              action_buttons1
            )
          } else {
            div(
              bslib::as_fill_carrier(),
              if(is.null(uploaded$contrasts.csv)){
                bslib::layout_columns(
                  col_widths = c(-3,6,-3),
                  ## row_heights = c(1,6,"auto"),
                  gap = "0.3rem",
                  div( bs_alert("The comparison file (comparisons.csv) is an optional input file. The file contains a list of pre-defined comparisons between groups (e.g. treatment versus controls, mutant versus wild-type). If you do not have a comparisons file, you can create comparisons using the interactive comparison builder.", closable = FALSE), style = "margin-bottom: -50px;"),
                  bslib::card(
                    fileInputArea(
                      ns("contrasts_csv"),
                      div(
                        style = "display: flex; flex-direction: column; gap: 10px;",
                        shiny::h4("Upload comparisons.csv (optional)", class='mb-0')
                      ),
                      multiple = FALSE,
                      accept = c(".csv"),
                      width = "100%"
                    ),
                    style = "background-color: aliceblue; border: 0.1rem dashed steelblue;"
                  ),
                  action_buttons2
                )
              } else {
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
                  ),
                  action_buttons2
                )
              }
            )
          }
        )
      }

    })

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
      
      # if contrasts not in file name, give warning and return
      if(!grepl("contrast", input$contrasts_csv$name, ignore.case = TRUE) && !grepl("comparison", input$contrasts_csv$name, ignore.case = TRUE)){
        shinyalert::shinyalert(
          title = "Comparison not in filename.",
          text = "Please make sure the file name contains 'comparison', such as comparison_dataset.csv or comparison.csv.",
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
