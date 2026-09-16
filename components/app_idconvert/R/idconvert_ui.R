##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

idconvert_ui <- function(id) {
  ns <- shiny::NS(id) ## namespace

  organism_choices <- convert_organism_choices()
  
  board <- bslib::layout_columns(
    col_widths = c(2, 10),
    bslib::layout_columns(
      col_widths = 12,
      fill = FALSE,
      gap = "0.7rem",
      shiny::selectizeInput(
        ns("organism"), "Organism:",
        choices = NULL,
        options = list(maxOptions = length(organism_choices))
      ),
      shiny::div(
        paste(length(organism_choices), "species available"),
        style = "font-size: 11px; color: #888; margin: -10px 0 0 0;"
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
          "font-size: 12px; color: #888; margin-top: -10px;",
          "margin-bottom: 12px; display: inline-block;",
          "text-decoration: underline;"
        ),
        shiny::actionLink(ns("example"), "Load example features",
          style = "margin-right: 15px;"),
        shiny::actionLink(ns("clear"), "Clear")
      ),
      shiny::checkboxInput(ns("human_ortholog"),"Human ortholog",TRUE),
      shiny::conditionalPanel(
        condition = "input.human_ortholog == false",
        ns = ns,
        shiny::selectInput(ns("ortholog"), "Ortholog species:", choices=NULL)
      ),
      br(),
      div(
        style = "display: flex; flex-direction: column; gap: 0;",
        shiny::actionButton(ns("convert"), "Convert",
          icon = icon("arrows-rotate"), class = "btn-primary mb-2",
          width = "100%"
        ),
        shiny::uiOutput(ns("download_ui"))
      )
    ),
    bslib::layout_columns(
      col_widths = 12,
      div(
        style = "padding-left: 30px; height: 100%;",
        shiny::uiOutput(ns("table_area"), style = "height: 100%;")
      )
    )
  )

  ui <- OmicsBoardUI(
    ns = ns,
    title = "ID Converter",
    board
  )

  return(ui)
}
