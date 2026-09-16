##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

idconvert_ui <- function(id) {
  ns <- shiny::NS(id) ## namespace

  organism_choices <- convert_organism_choices()

  ui <- bslib::layout_columns(
    col_widths = c(2, 10),
    class = "p-3",
    gap = "2rem",
    height = "calc(100vh - 72px)",
    bslib::layout_columns(
      col_widths = 12,
      gap = "14px",
      fill = FALSE,
      shiny::selectizeInput(ns("organism"), "Organism:",
        choices = c("Human" = "Human"),
        selected = "Human"
      ),
      shiny::div(
        paste(length(organism_choices), "species available"),
        style = "font-size: 12px; color: #888; margin: -12px 0 0 0;"
      ),
      shiny::selectInput(ns("datatype"), "Datatype:",
        choices = convert_datatype_choices(),
        selected = "proteomics"
      ),
      shiny::textAreaInput(ns("features"), "Gene/feature IDs (one per line):",
        rows = 20,
        placeholder = "e.g.\nTrp53\nENSMUSG00000059552\n..."
      ),
      div(
        style = paste(
          "font-size: 12px; color: #888; margin-top: -8px;",
          "margin-bottom: 0; display: inline-block;",
          "text-decoration: underline;"
        ),
        shiny::actionLink(ns("example"), "Load example features",
          style = "margin-right: 15px;"),
        shiny::actionLink(ns("clear"), "Clear")
      ),
      div(
        shiny::checkboxInput(ns("human_ortholog"),"Human ortholog",TRUE),
        shiny::conditionalPanel(
          condition = "input.human_ortholog == false",
          ns = ns,
          shiny::selectInput(ns("ortholog"), "Ortholog species:", choices=NULL)
        )
      ),
      br(),
      div(
        style = "display: flex; flex-direction: column; gap: 0;",
        shiny::actionButton(ns("convert"), "Convert",
          icon = icon("arrows-rotate"), class = "btn-primary mb-2",
          width = "100%"
        ),
        shiny::uiOutput(ns("download_ui"))
      ),
      br()
    ),
    bslib::layout_columns(
      col_widths = 12,
      div(
        style = "padding-left: 30px; height: 100%;",
        shiny::uiOutput(ns("table_area"), style = "height: 100%;")
      )
    )
  )

  title <- HTML("ID Converter <span style='font-size: 0.7em;'>&mdash; convert and annotate your features</span>")
  
  board <- OmicsBoardUI(
    ns = ns,
    #title = "AI Copilot",
    title = div(title, style="margin-left: 14px;"),
    info = FALSE,
    header_margin = "0px",
    ui
  )
  
  return(board)
}