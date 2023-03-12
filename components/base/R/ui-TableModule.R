tableWidget <- function(id) {
    ns <- shiny::NS(id)
    shiny::uiOutput(ns("widget"))
}

tableModule <- function(input, output, session,
                        func, func2=NULL, info.text="Info text",
                        title=NULL, label="", server=TRUE,
                        caption=NULL, caption2=caption,
                        csvFunc=NULL, filename="data.csv", ##inputs=NULL,
                        width=c("100%","100%"), height=c("auto","auto"),
                        options = NULL, info.width="300px",
                        selector = c("none","single", "multi","key")[1]
                        )
{
    ##require(bsutils)
    ns <- session$ns

    if(any(class(caption)=="reactive")) {
        caption <- caption()
    }
    if(class(caption)=="character") {
        caption <- shiny::HTML(caption)
    }

    options.button <- ""
    if(!is.null(options) && length(options)>0) {
        options.button <- shinyWidgets::dropdownButton(
            options,
            ##shiny::br(),
            ##dload,
            circle = TRUE, size = "xs", ## status = "danger",
            ## icon = shiny::icon("gear"),
            icon = shiny::icon("bars"),
            width = "250px",
            inputId = ns("options"),
            tooltip = shinyWidgets::tooltipOptions(title = "Settings", placement = "right")
        )
    }

    ##if(!is.null(label) && label!="") label <- paste0("(",label,")")
    label1 = shiny::HTML(paste0("<span class='module-label'>",label,"</span>"))
    title1 = title
    if(label!="") {
        title1 = shiny::HTML(paste0(title," (",label,")"))
    }

    zoom.button <- modalTrigger(
        ns("zoombutton"),
        ns("tablePopup"),
        icon("window-maximize"),
        class="btn-circle-xs"
    )

    header <- shiny::fillRow(
        ##flex=c(NA,NA,NA,NA,1),
        flex=c(NA,1,NA,NA,NA,NA),
        label1,
        shiny::div(class='plotmodule-title', title=title, title1),
        shinyWidgets::dropdownButton(
            shiny::tags$p(shiny::HTML(info.text)),
            shiny::br(),
            circle = TRUE, size = "xs", ## status = "danger",
            icon = shiny::icon("info"), width = info.width,
            inputId = ns("info"), right=FALSE,
            tooltip = shinyWidgets::tooltipOptions(title = "Info", placement = "right")
        ),
        options.button,
        shiny::div(class='download-button',
            shinyWidgets::dropdownButton(
                shiny::downloadButton(ns("csv"), "CSV"),
                circle = TRUE, size = "xs", ## status = "danger",
                icon = shiny::icon("download"), width = "80px", right=FALSE,
                tooltip = shinyWidgets::tooltipOptions(title = "Download",
                    placement = "right")
            )),
        ##withTooltip(zoom.button,"maximize table")
        zoom.button
    )

    CSVFILE = paste0(gsub("file","data",tempfile()),".csv")
    CSVFILE

    ## render2 <- shiny::renderPlot({plot_array[[3]]()}, res=res)
    download.csv <- shiny::downloadHandler(
        filename = filename,
        content = function(file) {
          if(!is.null(csvFunc)) {
            dt <- csvFunc()
          } else {
            dt <- func()$x$data
          }
          ##write.csv(dt, file=CSVFILE, row.names=FALSE)
          ##file.copy(CSVFILE, file, overwrite=TRUE)
          write.csv(dt, file=file, row.names=FALSE)
        }
    )
    output$csv <- download.csv

    if(is.null(func2)) func2 <- func
    if(length(height)==1) height <- c(height,height)
    if(length(width)==1)  width  <- c(width,width)
    ##ifnotchar.int <- function(s) ifelse(grepl("[%]$|auto|vmin|vh|vw|vmax",s),s,as.integer(s))
    ifnotchar.int <- function(s) suppressWarnings(
      ifelse(!is.na(as.integer(s)), paste0(as.integer(s),"px"), s))
    width.1  <- ifnotchar.int(width[1])
    width.2  <- ifnotchar.int(width[2])
    height.1 <- ifnotchar.int(height[1])
    height.2 <- ifnotchar.int(height[2])

    output$datatable <- DT::renderDT({
        # If the options `scrollX` or `autoWidth` or `selector` are set,
        # the global defaults of the global.R
        # will be overwritten. This ensures those options
        # are kept so that the header scrolls properly, and clickable
        # properties for tables.
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
        # Remove striping and borders from all tables
        dt$x$container <- stringr::str_remove(dt$x$container, "stripe")
        dt$x$container <- stringr::str_remove(dt$x$container, "table-bordered")
        dt
    },
    fillContainer = T)

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
        dt$x$container <- stringr::str_remove(dt$x$container, "stripe")
        dt$x$container <- stringr::str_remove(dt$x$container, "table-bordered")
        dt
    },
    fillContainer = T)

    output$popuptable <- shiny::renderUI({
        if(any(class(caption2)=="reactive")) {
            caption2 <- caption2()
        }
        if(any(class(caption2)=="character")) {
            caption2 <- shiny::HTML(caption2)
        }
        shiny::tagList(
            shiny::div( caption2, class="caption2"),
            DT::DTOutput(ns("datatable2"), width=width.2, height=height.2)
        )
    })

    output$widget <- shiny::renderUI({

        modaldialog.style <- paste0("#",ns("tablePopup")," .modal-dialog {width:",width.2,";}")
        modalbody.style   <- paste0("#",ns("tablePopup")," .modal-body {min-height:",height.2,";}")
        modalfooter.none  <- paste0("#",ns("tablePopup")," .modal-footer{display:none;}")
        div.caption <- NULL
        if(!is.null(caption)) div.caption <- shiny::div(caption, class="table-caption")

        div(class="tablewidget",
            shiny::fillCol(
              flex = c(NA,NA,1,NA),
              shiny::tags$head(shiny::tags$style(modaldialog.style)),
              shiny::tags$head(shiny::tags$style(modalbody.style)),
              shiny::tags$head(shiny::tags$style(modalfooter.none)),
              div(header, class="plotmodule-header"),
              div.caption,
              DT::DTOutput(ns("datatable"), width=width.1, height=height.1),
              shiny::div(class="popup-table",
                         modalUI(
                           id = ns("tablePopup"),
                           title = title,
                           size = "fullscreen",
                           shiny::uiOutput(ns("popuptable"))
                         ))
            ))
        })

    module <- list(
        ##data = func,
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
}
