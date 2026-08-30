##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


StudioUI <- function(id) {
  ns <- shiny::NS(id)

  sysprompt <- "You are computational biologist and will be asked questions about the following experiment. Make specific observations if you can."

  studio_card <- bslib::navset_underline(
    bslib::nav_panel(
      "Studio",
      shiny::br(),
      ui <- bslib::layout_columns(
        col_widths = c(6,6),
        row_heights = "26px",
        shiny::actionButton(ns("show_reports"), "Reports", width = "100%",
          class = "xbtn", icon = icon("file-lines")),
        shiny::actionButton(ns("show_infographics"), "Infographics", width = "100%",
          class = "xbtn", icon = icon("image")),
        shiny::actionButton(ns("show_poster"), "Poster", width = "100%",
          class = "xbtn", icon = icon("note-sticky")),
        shiny::actionButton(ns("show_slidedeck"), "Slide deck", width = "100%",
          class = "xbtn", icon = icon("film"))
      )
    )
  )

  settings_card <- bslib::navset_hidden(
    id = ns("settings"),
    selected = "Reports",
    bslib::nav_panel(
      "Reports",
      AiReportSettings(ns("aireport"))
    ),
    bslib::nav_panel(
      "Infographics",
      InfographicSettings(ns("infographic"))
    ),
    bslib::nav_panel(
      "Poster",
      VisReportSettings(ns("poster"), output_format="poster")
    ),    
    bslib::nav_panel(
      "Slide deck",
      VisReportSettings(ns("slide"), output_format="slide", type="comparison")
    )
  )

  output_panel <- bslib::navset_hidden(
    id = ns("studiopanel"),
    selected = "Reports",
    bslib::nav_panel(
      "Reports",
      AiReportUI(ns("aireport"))
    ),
    bslib::nav_panel(
      "Infographics",
      InfographicUI(ns("infographic"))
    ),
    bslib::nav_panel(
      "Poster",
      VisReportUI(ns("poster"), output_format="poster")
    ),    
    bslib::nav_panel(
      "Slide deck",
      VisReportUI(ns("slide"), output_format="slide", type="comparison")
    )
  )

  ui <- bslib::layout_columns(
    col_widths = c(3, 9),
    height = "calc(100vh - 76px)",    
    fill = TRUE,
    ## left sidebar
    bslib::layout_columns(
      col_widths = c(12),
      height = "100%",
      bslib::card(
        studio_card
      ),
      bslib::card(
        settings_card
      )
    ),
    bslib::card(
      output_panel
    )
  )

  board <- OmicsBoardUI(
    ns = ns,
    #title = HTML("AI Studio <small><small>&mdash; your creative space</small></small>"),
    title = HTML("Studio <span style='font-size: 0.7em;'>&mdash; your creative space</span>"),
    info = FALSE,
    header_margin = "0px",
    ui
  )
  
  return(board)
}
