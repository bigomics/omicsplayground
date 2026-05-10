##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


CopilotUI <- function(id) {
  ns <- shiny::NS(id)

  sysprompt <- "You are computational biologist and will be asked questions about the following experiment. Make specific observations if you can."

  examples_card <- bslib::navset_underline(
    bslib::nav_panel(
      "Examples",
      shiny::br(),
      "Example questions:",
      shiny::br(),
      shiny::actionButton(ns("ask_describe"), "Describe my experiment", width = "100%", class = "xbtn"),
      shiny::actionButton(ns("ask_findings"), "Summarize main findings", width = "100%", class = "xbtn"),
      shiny::actionButton(ns("ask_pathways"), "What pathways are involved?", width = "100%", class = "xbtn"),
      shiny::actionButton(ns("show_biomarkers"), "Show top biomarkers", width = "100%", class = "xbtn"),
      shiny::actionButton(ns("find_references"), "Find references", width = "100%", class = "xbtn"),
      shiny::actionButton(ns("get_expression"), "Get expression of MTOR", width = "100%", class = "xbtn"),
      shiny::actionButton(ns("plot_volcano"), "Show volcano plot", width = "100%", class = "xbtn"),
      shiny::br(),
      shiny::br(),
      shiny::checkboxInput(ns("followup"), "Suggest follow-up questions", TRUE)
    ),
    bslib::nav_panel(
      "Settings",
      br(),
      shiny::textAreaInput(ns("sysprompt"), "System prompt:", value = sysprompt, height = 120, width = "100%"),
      br(),
      shiny::radioButtons(ns("response_length"), "Response length:",
        choices = c("default","short","longer"), selected="short", inline = TRUE
      ),
      br(),
      shiny::checkboxGroupInput(ns("context"), "Context:",
        choices = NULL, inline = TRUE
      ),
      br(),
      actionButton(ns("reset"), "Reset model")
    )
  )
  
  chat_card <- bslib::card(
    class = "border-0",
    fill = FALSE,
    height = "calc(100vh - 100px)",
    bs_alert(HTML("<b>EXPERIMENTAL</b>. This AI module is experimental in early beta. Only for testing purposes."), closable = FALSE, style="warning"),
    shinychat::chat_ui(ns("chat"),
      width = "100%", height = "min(100%,770px)",
      fill = TRUE
    )
  )

  studio_card <- bslib::navset_underline(
    bslib::nav_panel(
      "Studio",
      shiny::br(),
      ui <- bslib::layout_columns(
        col_widths = c(6,6),
        row_heights = "26px",
        shiny::actionButton(ns("studio_podcast"), "Podcast", width = "100%",
          class = "xbtn", icon = icon("headphones")),
        shiny::actionButton(ns("studio_slidedeck"), "Slide deck", width = "100%",
          class = "xbtn", icon = icon("film")),
        shiny::actionButton(ns("studio_report"), "Reports", width = "100%",
          class = "xbtn", icon = icon("film")),
        shiny::actionButton(ns("studio_quiz"), "Quiz", width = "100%", class = "xbtn",
          icon = icon("quora")),
        shiny::actionButton(ns("studio_infographic"), "Infographic", width = "100%",
          class = "xbtn", icon = icon("image")),
        shiny::actionButton(ns("studio_datatable"), "Data table", width = "100%",
          class = "xbtn", icon = icon("table"))
      )
    )
  )

  output_card <- bslib::navset_underline(
    bslib::nav_panel(
      "Output",
      height = "600px",
      plotOutput(ns("plot"))
    )
  )

  ui <- bslib::layout_columns(
    col_widths = c(2, 7, 3),
    style = "height: min(90%,700)",
    fill = TRUE,
    ## left sidebar
    bslib::card(
      class = "border-0",
      fill = TRUE,
      height = "calc(100vh - 100px)",
      bslib::layout_columns(
        col_widths = c(12),
        examples_card      
      )
    ),
    ## center section
    chat_card,
    ## right sidebar
    bslib::layout_columns(
      col_widths = c(12),
      height = "calc(100vh - 100px)",
      bslib::card(
        class = "border-0",
        studio_card
      ),
      bslib::card(
        output_card
      )
    )
  )

  return(ui)
}

