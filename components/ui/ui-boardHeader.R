##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


header_infotext <- "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."

boardHeader <- function(title, info_link) {
  div(
    fillRow(
      flex = c(NA, 1),
      shiny::div(
        id = "navheader-current-section",
        HTML(paste0(title, " &nbsp;")),
        withTooltip(
          shiny::actionLink(
            inputId = info_link,
            label = "",
            icon = shiny::icon("youtube"),
            style = "color: #ccc;"
          ),
          "Video tutorial about this board"
        )
      ),
      div()
    )
  )
}

OmicsBoardUI <- function(id, title, ...) {
  ns <- shiny::NS(id) ## namespace
  bslib::layout_columns(
    col_widths = 12,
    gap = 0,
    row_heights = list("auto", 1),
    fillRow(
      flex = c(NA, 1, NA),
      style = "margin-top: -6px;",
      shiny::div(
        id = "navheader-current-section",
        HTML(paste0(title, " &nbsp;")),
        withTooltip(
          shiny::actionLink(
            inputId = ns("board_info"),
            label = "",
            icon = shiny::icon("youtube"),
            style = "color: #ccc;"
          ),
          "Video tutorial about this board"
        )
      ),
      div(),
      shiny::uiOutput(ns("current_dataset"))
    ),
    ...
  )
}

OmicsBoard <- function(id, pgx, title, infotext) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    output$current_dataset <- shiny::renderUI({
      has.pgx <- !is.null(pgx$name) && length(pgx$name) > 0
      shiny::req(has.pgx)
      pgx.name <- gsub(".*\\/|[.]pgx$", "", pgx$name)
      div(
        shiny::actionButton(
          ns("dataset_click"), pgx.name,
          class = "quick-button",
          style = "border: none; color: black; font-size: 1.2em; background-color: transparent;"
        ),
        style = "padding: 30px 0px 0px 0px;"
      )
    })
    
    ## ------- observe functions -----------
    shiny::observeEvent(input$board_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML(paste0("<strong>",title,"Board</strong>")),
        shiny::HTML(infotext),
        easyClose = TRUE, size = "xl"
      ))
    })

    ## Show experiment info if dataset name is clicked.
    observeEvent(input$dataset_click, {
      shiny::req(pgx$name)
      has.infographic <- !is.null(pgx$wgcna$report$infographic)
      pgx.name <- gsub(".*\\/|[.]pgx$", "", pgx$name)    
      if(has.infographic) {
        img <- pgx$wgcna$report$infographic
        footer <- gsub("- |\n"," ",pgx$wgcna$report$bullets)
        footer <- paste("<b>WGCNA graphical abstract</b>. ",footer)
        ui.showImageModal(img, title=NULL, footer, width=1088)         
      } else {
        fields <- c("name", "description", "datatype", "date", "settings",
          "omicsplayground_version")
        body <- playbase::pgx.info(pgx, fields=fields, format="html")
        shiny::showModal(shiny::modalDialog(
          header = pgx.name,
          div(HTML(body), style = "font-size: 1.1em;"),
          footer = NULL,
          size = "l",
          easyClose = TRUE,
          fade = FALSE
        ))
      }
    })
    
    
  })
}

