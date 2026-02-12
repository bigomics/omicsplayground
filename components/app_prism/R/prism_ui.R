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
    div(class = "navbar navbar-static-top", div(title, class = "container-fluid")),
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
          "minimal","classic","xkcd","prism")), selected="gray"),
        bslib::layout_columns(
          col_widths = c(6,6),
          shiny::sliderInput(ns("pointsize"), "Point size:", 1, 10, 3, step=1),
          shiny::sliderInput(ns("fontsize"), "Font size:", 8, 32, 18, step=2)
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
            shiny::plotOutput(ns("plot1"), height='600px'),
            ## shinychat::chat_ui(
            ##   ns("chartbot"),
            ##   style = "max-height: 180px; width: min(800px, 100%);",
            ##   fill = FALSE
            ## )
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
