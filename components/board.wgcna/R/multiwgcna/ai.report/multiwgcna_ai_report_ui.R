##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

multiwgcna_ai_report_ui <- function(id) {
  ns <- shiny::NS(id)
  controls_ns <- shiny::NS(ns("controls"))

  show_prompt_input <- shiny::checkboxInput(
    controls_ns("show_prompt"),
    "Show prompt",
    FALSE
  )

  infographic_options <- shiny::tagList(
    shiny::tags$hr(),
    shiny::selectInput(
      controls_ns("image_style"),
      "Infographic Style:",
      choices = NULL,
      width = "100%"
    ),
    shiny::radioButtons(
      controls_ns("image_blocks"),
      "Layout:",
      choices = c("1 Panel" = "1", "2 Panels" = "2", "3 Panels" = "3"),
      selected = "1",
      inline = TRUE
    )
  )

  multiomics_ai_report_layout_ui(
    ns("layout"),
    text_title = "AI MultiWGCNA Report",
    diagram_title = "Cross-Layer Module Diagram",
    infographic_title = "Graphical Abstract",
    text_options = show_prompt_input,
    infographic_options = infographic_options
  )
}

multiwgcna_ai_report_inputs_ui <- function(id) {
  ns <- shiny::NS(id)
  multiwgcna_ai_report_controls_ui(ns("controls"))
}

multiwgcna_ai_report_layout_server <- function(id, text_reactive, watermark = FALSE) {
  multiomics_ai_report_layout_server(
    id,
    text_reactive = text_reactive,
    filename = "multiwgcna-ai-report",
    title = "MultiWGCNA AI Report",
    watermark = watermark
  )
}
