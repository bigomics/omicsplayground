##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

convert_ui <- function(id) {
  ns <- shiny::NS(id) ## namespace

  title <- div("Gene ID Converter", style = "font-size: 18px;")

  organism_choices <- convert_organism_choices()

  ui <- bslib::page_fillable(
    padding = 0,
    div(class = "navbar navbar-static-top", div(title, class = "container-fluid"),
      style = "margin-top: 24px;"),
    bslib::layout_columns(
      col_widths = c(2, 10),
      class = "p-3",
      gap = "2rem",
      height = "calc(100vh - 72px)",
      bslib::layout_columns(
        col_widths = 12,
        fill = FALSE,
        shiny::selectizeInput(ns("organism"), "Organism:",
          choices = c("Human" = "Human"),
          selected = "Human"
        ),
        shiny::div(
          paste(length(organism_choices), "species available"),
          style = "font-size: 11px; color: #888; margin-top: -12px;"
        ),
        shiny::selectInput(ns("datatype"), "Datatype:",
          choices = convert_datatype_choices(),
          selected = "proteomics"
        ),
        shiny::textAreaInput(ns("features"), "Gene/feature IDs (one per line):",
          rows = 20,
          placeholder = "e.g.\nTrp53\nENSMUSG00000059552\n..."
        ),
        shiny::actionLink(ns("example"), "Load example features",
          style = paste(
            "font-size: 12px; color: #888; margin-top: -12px;",
            "margin-bottom: 12px; display: inline-block;",
            "text-decoration: underline;"
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
  )

  return(ui)
}
