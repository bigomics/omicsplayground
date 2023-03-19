##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics Sagl. All rights reserved.
##


TableModuleUI <- function(id,
                          height = c(400,800),
                          width = c("auto","100%"),
                          info.text="Figure",
                          title="",
                          options = NULL,
                          label="",
                          caption="",
                          caption2=info.text,
                          just.info=FALSE,
                          show.maximize = TRUE) {
  ns <- shiny::NS(id)

  if(length(height)==1) height <- c(height,800)
  if(length(width)==1)  width  <- c(width,1200)

  ifnotchar.int <- function(s) suppressWarnings(
    ifelse(!is.na(as.integer(s)), paste0(as.integer(s),"px"), s))
  width.1  <- ifnotchar.int(width[1])
  width.2  <- ifnotchar.int(width[2])
  height.1 <- ifnotchar.int(height[1])
  height.2 <- ifnotchar.int(height[2])

  options.button <- ""

  if(!just.info && !is.null(options) && length(options)>0) {
    options.button <- DropdownMenu(
      options,
      size = "xs",
      icon = shiny::icon("bars"),
      status = "default"
    )
  }

  dload.button <- DropdownMenu(
    div(
      style = "width: 150px;",
      shiny::a("Download table data",
               style = "text-align: center;"),
      shiny::br(),
      shiny::hr(),
      div(style = "text-align: center;",
        shiny::downloadButton(
          outputId = ns("download"),
          label = "Download",
        )
      )
    ),
    size = "xs",
    icon = shiny::icon("download"),
    status = "default"
  )

  if(!is.null(label) && label!="") label <- paste0("&nbsp;(",label,")")
  label <- shiny::div(class = "tablemodule-title", shiny::HTML(label))

  zoom.button <- NULL
  if(show.maximize) {
    zoom.button <- modalTrigger(ns("zoombutton"),
                                ns("datatablePopup"),
                                icon("window-maximize"),
                                class="btn-circle-xs"
    )
  }

  header <- shiny::fillRow(
    flex = c(NA,1,NA,NA,NA,NA),
    shiny::div(class='tablemodule-title', title=title, title),
    label,
    DropdownMenu(
      shiny::div(class='tablemodule-info', shiny::HTML(info.text)),
      width = "250px",      
      size = "xs",
      icon = shiny::icon("info"),
      status = "default"
    ),
    options.button,
    shiny::div(class='download-button', title='download', dload.button),
    shiny::div(class='zoom-button', title='zoom', zoom.button)
  )

  # Modal stuff

  popupdatatableUI <- function() {
    if(any(class(caption2)=="reactive")) {
      caption2 <- caption2()
    }
    if(any(class(caption2)=="character")) {
      caption2 <- shiny::HTML(caption2)
      caption2 <- shiny::div(caption2, class="caption2 popup-table-caption")
    }
    shiny::tagList(
      shiny::div(
        class = "popup-table",
        DT::DTOutput(ns("datatable2"), width=width.2, height=height.2)
      ),
      caption2
    )
  }
  modaldialog.style <- paste0("#",ns("plotPopup")," .modal-dialog {width:",width.2,";}")
  modalbody.style <- paste0("#",ns("plotPopup")," .modal-body {min-height:",height.2,"; padding:30px 300px;}")
  modalcontent.style <- paste0("#",ns("plotPopup")," .modal-content {width:100vw;}")
  modalfooter.none <- paste0("#",ns("plotPopup")," .modal-footer{display:none;}")

  # Div construction
  div( class="tablemodule",
       shiny::fillCol(
         flex = c(NA,1,NA,0.001,NA),
         height = height.1,
         div( header, class="tablemodule-header"),
         DT::DTOutput(ns("datatable"), width=width.1, height=height.1) %>%
           shinycssloaders::withSpinner(),
         div(
           class = "footer",
           shiny::HTML(caption)
         ),
         shiny::div(class="popup-modal",
                    modalUI(
                      id = ns("datatablePopup"),
                      title = div(class="popup-table-title", title),
                      size = "fullscreen",
                      footer = NULL,
                      popupdatatableUI()
                    )
         ),
         shiny::tagList(
           shiny::tags$head(shiny::tags$style(modaldialog.style)),
           shiny::tags$head(shiny::tags$style(modalbody.style)),
           shiny::tags$head(shiny::tags$style(modalcontent.style)),
           shiny::tags$head(shiny::tags$style(modalfooter.none))
         )
       )
  )
}

TableModuleServer <- function(id,
                              func,
                              func2 = NULL,
                              csvFunc = NULL,
                              height = c(640,800),
                              width = c("auto",1400),
                              selector = c("none","single", "multi","key")[1],
                              filename = "data.csv")
{
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      if(is.null(func2)) func2 <- func

      # Downloader
      output$download <- shiny::downloadHandler(
        filename = filename,
        content = function(file) {
          if(!is.null(csvFunc)) {
            dt <- csvFunc()
          } else {
            dt <- func()$x$data
          }
          dt2 <- format(dt, digits=4)  ## round
          write.csv(dt2, file=file, row.names=FALSE)
        }
      )

      output$datatable <- DT::renderDT({
        # If the options `scrollX` or `autoWidth`, `fillContainer` or
        # `selector` are set, the global defaults of the global.R will
        # be overwritten. This ensures those options are kept so that
        # the header scrolls properly, and clickable properties for
        # tables.
        dt <- func()
        active_options <- names(dt$x$options)
        if("scrollX" %in% active_options){
          dt$x$options$scrollX <- TRUE
        }
        if("autoWidth" %in% active_options){
          dt$x$options$autoWidth <- FALSE
        }
        if(!is.null(selector)){
          dt$x$selection$mode = selector
        }
        if(!is.null(dt$x$fillContainer)){
          dt$x$fillContainer = FALSE
        }
        # Remove striping and borders from all tables
        dt$x$container <- stringr::str_remove(dt$x$container, "stripe")
        dt$x$container <- stringr::str_remove(dt$x$container, "table-bordered")
        dt
      },
      fillContainer = TRUE)

      output$datatable2 <- DT::renderDT({
        dt <- func2()
        active_options <- names(dt$x$options)
        if("scrollX" %in% active_options){
          dt$x$options$scrollX <- TRUE
        }
        if("autoWidth" %in% active_options){
          dt$x$options$autoWidth <- FALSE
        }
        if(!is.null(selector)){
          dt$x$selection$mode = selector
        }
        if(!is.null(dt$x$fillContainer)){
          dt$x$fillContainer = FALSE
        }
        dt$x$container <- stringr::str_remove(dt$x$container, "stripe")
        dt$x$container <- stringr::str_remove(dt$x$container, "table-bordered")
        dt
      },
      fillContainer = TRUE)

      module <- list(
        data = shiny::reactive(func()$x$data),
        rows_current = shiny::reactive(input$datatable_rows_current),
        rows_selected = shiny::reactive(input$datatable_rows_selected),
        rows_all = shiny::reactive(input$datatable_rows_all),
        rownames_current = shiny::reactive({
          rns <- rownames(func()$x$data)
          if(is.null(rns)) rns <- 1:nrow(func()$x$data)
          rns[input$datatable_rows_current]
        }),
        rownames_selected = shiny::reactive({
          rns <- rownames(func()$x$data)
          if(is.null(rns)) rns <- 1:nrow(func()$x$data)
          rns[input$datatable_rows_selected]
        }),
        rownames_all = shiny::reactive({
          rns <- rownames(func()$x$data)
          if(is.null(rns)) rns <- 1:nrow(func()$x$data)
          rns[input$datatable_rows_all]
        })
      )
      return(module)
    })
}
