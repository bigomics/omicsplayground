##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## used for navset_pill in main. fullheight-page is on by default so every
## panel fills to the window bottom instead of sizing to content -- pass
## class = NULL to opt out.
omicspanel <- function(p, class = "fullheight-page") {
  div(p, style="margin: 0 12px 0 12px;", class = trimws(paste("omicspanel", class)))
}

## header_margin: margin-top of the title row. Boards that need extra
## breathing room at the top of the window (AI Studio, Obi AI) pass a
## positive value; everything else keeps the -12px that tucks the title
## baseline in line with the sidebar's "Menu" heading.
OmicsBoardUI <- function(ns, title, ..., subtitle=NULL, info=TRUE, header_margin="-12px") {
  div.info <- NULL
  if(isTRUE(info)) {
    div.info <- withTooltip(
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
    height = "100%",
    class = "omicsboard",
    style = "padding-bottom: 5px;",
    row_heights = list("auto", "auto", 1),
    fillRow(
      flex = c(NA, 1, NA),
      style = paste0("margin-top: ", header_margin, ";"),
      shiny::div(
        id = "navheader-current-section",
        HTML(paste0(title, "&nbsp;")),
        div.info
      ),
      shiny::uiOutput(ns("current_dataset")),
      div(shiny::actionButton(ns("logo_click"),
        shiny::tags$img(src = "static/bigomics-logo-small.png", height = "28px"),
        class = "quick-button", style="border: 0px; background-color: transparent;"),
        style="margin: 6px 5px -20px 0; padding: 0 0 0 80px;")
    ),
    div(bigdash::bd_visibility_probe(ns), height="0px", style='width: 100%;'),
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
      ##pgx.name <- "(no dataset)"
      pgx.name <- ""
    }
    div(
      shiny::actionButton(
        ns("dataset_click"), pgx.name,
        class = "quick-button",
        style = "border: none; color: black; font-size: 1.2em; background-color: transparent;"
      ),
      style = "padding: 13px 0px 0px 0px; text-align: center;"
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
  ##
  ## The infographic now lives under the v1 AI-report schema at
  ## pgx$ai$<slot>$infographic (playbase::pgx.update_infographics(),
  ## ai-infographics.R) -- "combined" is the dataset-wide report, labelled
  ## "Summary" in AI Studio's own infographic tab (.AI_INFOGRAPHIC_LABELS).
  ## $infographic holds raw PNG bytes + a status, not the decoded raster
  ## array ui.showImageModal() expects, so decode via png::readPNG() first.
  ## The old schema's short "bullets" field (a few summary lines meant for
  ## the modal footer) has no equivalent here -- pgx$ai$combined$report is a
  ## full multi-section markdown write-up, not caption-sized -- so the
  ## footer is just a static label now.
  observeEvent(input$dataset_click, {
    shiny::req(pgx$name)
    pgx.name <- gsub(".*\\/|[.]pgx$", "", pgx$name)
    infographic <- pgx$ai$combined$infographic
    has.infographic <- !is.null(infographic) &&
      identical(infographic$status, "done") &&
      length(infographic$bytes) > 0
    img <- if (has.infographic) {
      tryCatch(png::readPNG(infographic$bytes), error = function(e) NULL)
    } else {
      NULL
    }
    if(!is.null(img)) {
      ui.showImageModal(img, title=NULL, "<b>Summary infographic</b>", width=1088)
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

