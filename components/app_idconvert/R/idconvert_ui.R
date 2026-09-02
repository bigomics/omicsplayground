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
      div(
        style = "margin: -15px 0 0px 0;",             
        shiny::textAreaInput(ns("features"), "Gene/feature IDs (one per line):",
          rows = 16,
          placeholder = "e.g.\nTrp53\nENSMUSG00000059552\n..."
        ),
        shiny::actionLink(ns("example"), "Load example features",
          style = paste(
            "font-size: 12px; color: #888; margin-top: -15px;",
            "margin-bottom: 8px; display: inline-block;",
            "text-decoration: underline; padding-right: 10px;"
          )
        ),
        shiny::actionLink(ns("clear"), "Clear",
          style = paste(
            "font-size: 12px; color: #888; margin-top: -12px;",
            "margin-bottom: 8px; display: inline-block;",
            "text-decoration: underline;"
          )
        )
      ),
      div(
        style = "margin: -5px 0 10px 0;",     
        shiny::checkboxInput( ns("use_bridge"), "Use bridge organism", FALSE),
        shiny::conditionalPanel(
          condition = "input.use_bridge",
          ns = ns,
          shiny::selectizeInput(ns("bridge_organism"), "Bridge organism(s):",
            choices = NULL, multiple = TRUE
          )
        )
      ),
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
