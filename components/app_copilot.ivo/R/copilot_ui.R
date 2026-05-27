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
      shiny::br(),
      shiny::actionButton(ns("ask_findings"), "Summarize main findings", width = "100%", class = "xbtn"),
      shiny::actionButton(ns("ask_pathways"), "What pathways are involved?", width = "100%", class = "xbtn"),
      shiny::actionButton(ns("show_biomarkers"), "Show top biomarkers", width = "100%", class = "xbtn"),
      shiny::actionButton(ns("find_references"), "Find references", width = "100%", class = "xbtn"),
      shiny::actionButton(ns("get_expression"), "Get expression of MTOR", width = "100%", class = "xbtn")
    ),
    bslib::nav_panel(
      "Settings",
      br(),
      shiny::textAreaInput(ns("sysprompt"), "System prompt:", value = sysprompt, height = 100, width = "100%"),
      br(),
      shiny::radioButtons(ns("response_length"), "Response length:",
        choices = c("shorter","longer"), selected="shorter", inline = TRUE
      ),
      br(),      
      shiny::checkboxInput(ns("followup"), "Suggest follow-up questions", TRUE),
      br(),      
      actionButton(ns("reset"), "Reset model")
    )
  )

  input_card <- bslib::navset_underline(
    bslib::nav_panel(
      "Input sources",
      br(),
      shiny::checkboxGroupInput(ns("context"), "Dataset context:",
        choices = NULL, inline = TRUE
      ),
      br(),
      br(),      
      shiny::fileInput("file", "Add sources:", accept = NULL, multiple=TRUE),
    )
  )

  chat_card <- bslib::card(
    class = "border-0",
    fill = FALSE,
    height = "calc(100vh - 80px)",
    bs_alert(HTML("<b>EXPERIMENTAL</b>. This AI module is experimental in early beta. Only for testing purposes."), closable = FALSE, style="warning"),
    shinychat::chat_ui(ns("chat"),
      width = "100%", height = "min(100%,770px)",
      fill = TRUE
    )
  )

  output_card <- bslib::navset_underline(
    bslib::nav_panel(
      "Output",
      height = "600px",
      plotOutput(ns("plot"))
    )
  )

  info_card <- bslib::navset_underline(
    bslib::nav_panel(
      "Info",
      height = "600px",
      div(style="height: 300px;")
    )
  )

  ui <- bslib::layout_columns(
    col_widths = c(3, 6, 3),
    style = "height: min(90%,700)",
    fill = TRUE,
    ## left sidebar
    bslib::layout_columns(
      col_widths = c(12),
      height = "calc(100vh - 80px)",
      bslib::card(
        examples_card
      ),
      bslib::card(
        input_card
      )
    ),
    ## center section
    chat_card,
    ## right sidebar
    bslib::layout_columns(
      col_widths = c(12),
      height = "calc(100vh - 80px)",
      bslib::card(
        info_card
      ),
      bslib::card(
        output_card
      )
    )
  )

  board <- OmicsBoardUI(
    id = ns("board"),
    #title = "AI Copilot",
    title = "ObiOne Copilot",    
    ui
  )
  
  return(board)
}

