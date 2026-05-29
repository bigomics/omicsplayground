##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


CopilotUI <- function(id) {
  ns <- shiny::NS(id)

  SYSPROMPT <- "You are computational biologist and will be asked questions about the following omics experiment."

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
    )
  )
  
  input_card <- bslib::navset_underline(
    bslib::nav_panel(
      "Settings",
      br(),
      shiny::selectInput(ns("style"), "Style:",
        choices = c("scientist","teacher","poet"), selected="scientist", 
      ),
      shiny::textAreaInput(ns("sysprompt"), NULL, value = "",
        height = 80, width = "100%"),
      shiny::radioButtons(ns("response_length"), "Response length:",
        choices = c("shorter","longer"), selected="shorter", inline = TRUE
      ),
      shiny::checkboxInput(ns("followup"), "Suggest follow-up questions", TRUE),
      br(),      
      actionButton(ns("reset"), "Apply")
    ),
    bslib::nav_panel(
      "Context",
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
      ##div(style="height: 300px;")
      copilot_info_ui(ns("info"))
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
    title = "Obi-One Copilot",    
    ui
  )
  
  return(board)
}

