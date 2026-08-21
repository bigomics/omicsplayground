##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## used for navset_pill in main
omicspanel <- function(p) {
  div(p, style="margin: 0 8px 0 8px;", class = "omicspanel")
}  

header_infotext <- "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."

boardHeader <- function(title, info_link, mid=NULL, right=NULL) {
  if(is.null(mid)) mid <- div()
  if(is.null(right)) right <- div()  
  div(
    fillRow(
      flex = c(NA, 1, NA),
      shiny::div(
        id = "navheader-current-section",
        HTML(paste0(title, "&nbsp;")),
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
      mid,
      right
    )
  )
}

OmicsBoardUI <- function(ns, title, ..., info=TRUE) {
  div.link <- NULL
  if(isTRUE(info)) {
    div.link <- withTooltip(
      shiny::actionLink(
        inputId = ns("board_info"),
        label = "",
        icon = shiny::icon("youtube"),
        style = "color: #ccc;"
      ),
      "Video tutorial about this board"
    )
  }

  bslib::layout_columns(
    col_widths = 12,
    gap = 0,
    row_heights = list("auto", 1),
    fillRow(
      flex = c(NA, 1, NA),
      style = "margin-top: -12px;",
      shiny::div(
        id = "navheader-current-section",
        HTML(paste0(title, "&nbsp;")),
        div.link
      ),
      shiny::uiOutput(ns("current_dataset")),
#      div( tags$a(
#        href="http://www.bigomics.ch", target = "_blank", 
#        shiny::tags$img(src = "static/bigomics-logo-small.png", height = "28px")),
#        style="margin: 34px 5px -20px 0; padding: 0 0 0 80px;")
      div(shiny::actionButton(ns("logo_click"),
        shiny::tags$img(src = "static/bigomics-logo-small.png", height = "28px"),
        class = "quick-button", style="border: 0px; background-color: transparent;"),
        style="margin: 28px 5px -20px 0; padding: 0 0 0 80px;")
    ),
    bigdash::bd_visibility_probe(ns),
    div(...)
  )
}

OmicsBoard <- function(session, pgx, title, infotext=NULL, purge=NULL) {
  input <- session$input
  output <- session$output
  ns <- session$ns ## NAMESPACE

  bigdash::bd_is_visible(
    input,
    purge = if (is.null(purge)) isTRUE(opt$ENABLE_PLOTLY_PURGE) else purge
  )

  output$current_dataset <- shiny::renderUI({
    has.pgx <- !is.null(pgx$name) && length(pgx$name) > 0
    if(has.pgx) {
      pgx.name <- gsub(".*\\/|[.]pgx$", "", pgx$name)
    } else {
      pgx.name <- "(no dataset)"
    }
    div(
      shiny::actionButton(
        ns("dataset_click"), pgx.name,
        class = "quick-button",
        style = "border: none; color: black; font-size: 1.2em; background-color: transparent;"
      ),
      style = "padding: 30px 0px 0px 0px; text-align: center;"
    )
  })

  ## ------- observe functions -----------
  shiny::observeEvent(input$board_info, {
    if(!is.null(infotext)) {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML(paste0("<strong>",title,"&nbsp;Board</strong>")),
        shiny::HTML(infotext),
        easyClose = TRUE, size = "xl"
      ))
    }
  })

  shiny::observeEvent(input$logo_click, {
    ui.showAboutModal()
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
}

