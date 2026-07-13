##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

prism_ui <- function(id) {
  ns <- shiny::NS(id) ## namespace

  buttons <- div(
    class = "toolbox",
    bslib::layout_columns(
      col_widths = c(3,3,3,3),
      gap = 0,
      actionButton(ns("button1"), "button1", icon=icon("pen")),
      actionButton(ns("button2"), "button2", icon=icon("eraser")),
      actionButton(ns("button3"), "button3", icon=icon("font")),
      actionButton(ns("button4"), "button4", icon=icon("eye-dropper"))
    )
  )

  title <- div("BigOmics Prism", style="font-size: 18px;")
  
  ui <- bslib::page_fluid(
    tags$head(
      tags$script(type = "module", src = "static/prism-webr.js"),
      tags$style(HTML("
        .prism-chat-messages {
          max-height: 160px; overflow-y: auto; padding: 8px;
          display: flex; flex-direction: column; gap: 6px;
          border: 1px solid #e2e8f0; border-radius: 6px; margin-bottom: 8px;
        }
        .prism-msg {
          padding: 8px 12px; border-radius: 10px; max-width: 90%;
          font-size: 0.85rem; line-height: 1.4;
          white-space: pre-wrap; word-break: break-word;
        }
        .prism-msg-user { background:#6366f1; color:#fff; align-self:flex-end; }
        .prism-msg-assistant { background:#f1f5f9; color:#334155; align-self:flex-start; border:1px solid #e2e8f0; }
        .prism-msg-error { background:#fee2e2; color:#991b1b; align-self:flex-start; border:1px solid #fca5a5; }
        .prism-plot-error {
          color:#dc2626; background:#fee2e2; padding:1rem; border-radius:8px;
          font-family:monospace; font-size:0.82rem; display:none; margin-top:8px;
        }
        .prism-webr-status {
          font-size:0.8rem; padding:5px 10px; border-radius:5px;
          background:#f8fafc; border:1px solid #e2e8f0; color:#64748b; margin-bottom:8px;
        }
      "))
    ),
    div(class = "navbar navbar-static-top", div(title, class = "container-fluid"),
      style="margin-top: 24px;"),
    bslib::layout_columns(
      col_widths = c(3,9),
      class = "p-3",
      bslib::layout_columns(
        col_widths = 12,
        fill = FALSE,
        buttons,
        shiny::selectInput(ns("dataset"), "Dataset:", c("mtcars","iris","geiger"),
          selected="mtcars"),
        shiny::selectInput(ns("theme"), "Theme:", sort(c("gray","bw","light","dark",
          "minimal","classic","prism")), selected="gray"),
        bslib::layout_columns(
          col_widths = c(6,6),
          shiny::sliderInput(ns("pointsize"), "Point size:", 1, 10, 3, step=1),
          shiny::sliderInput(ns("fontsize"), "Font size:", 8, 48, 18, step=4)
        ),
        wellPanel(
          style = "width: 100%; font-family: monospace; font-size: 11px;",
          shiny::htmlOutput(ns("plotcode"), height="400px")
        )
      ),
      bslib::layout_columns(
        col_widths = 12,
        class = "pl-4",
        bslib::navset_tab(
          bslib::nav_panel(
            title = "plot",
            div(id = ns("webr-status"), class = "prism-webr-status",
              "Initializing webR runtime…"),
            div(id = ns("plot-placeholder"),
              style = "color:#94a3b8; text-align:center; padding:2rem;",
              tags$div(style = "font-size:2.5rem;", "\U0001f4c8"),
              tags$div("Your plot will appear here")
            ),
            div(id = ns("plot-container")),
            div(id = ns("plot-error"), class = "prism-plot-error"),
            div(id = ns("chat-messages"), class = "prism-chat-messages"),
            div(
              shiny::textInput(ns("chartbot_user_input"),"",
                placeholder = "What do you want to plot?", width=600),
              shiny::actionButton(ns("chartbot_send"),"send",
                icon = icon("arrow-right-from-bracket"))
            )
          ),
          bslib::nav_panel(
            title = "data",
            shiny::dataTableOutput(ns("data1"))
          )
        )
      )
    )
  )
  
  return(ui)
}
