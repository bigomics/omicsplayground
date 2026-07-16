##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' Split a pasted block of gene/feature IDs into a character vector.
#' Splits on newlines, commas and other whitespace; trims and drops blanks.
parse_feature_list <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) {
    return(character(0))
  }
  parts <- strsplit(text, "[,\\s]+", perl = TRUE)[[1]]
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

#' Build the organism choices for the selector from playbase's species table.
#' Names are the display labels shown in the dropdown; values are the
#' species names playbase::getProbeAnnotation() expects as `organism`.
convert_organism_choices <- function() {
  species_name <- playbase::allSpecies(col = "species_name")
  display_name <- playbase::allSpecies(col = "display_name")
  keep <- !is.na(species_name) & species_name != "No organism"
  species_name <- species_name[keep]
  display_name <- display_name[keep]
  names(species_name) <- display_name
  species_name[order(display_name)]
}

#' Datatype choices for the selector; values are what
#' playbase::getProbeAnnotation() expects as `datatype`.
convert_datatype_choices <- function() {
  c(
    "RNA-seq" = "RNA-seq",
    "Proteomics" = "proteomics",
    "Metabolomics" = "metabolomics",
    "Lipidomics" = "lipidomics",
    "Methylomics" = "methylomics",
    "Multi-omics" = "multi-omics"
  )
}

#' The application server-side logic
#'
#' @param id Shiny module id.
#' @export
convert_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {

    shiny::updateSelectizeInput(session, "organism",
      choices = convert_organism_choices(), selected = "Human", server = TRUE
    )

    shiny::observeEvent(input$example, {
      features <- shiny::withProgress(
        message = "Loading example genes...",
        playbase::getExampleFeatures(organism = input$organism, n = 20)
      )
      if (length(features) == 0) {
        shiny::showNotification(
          "Could not find example genes for the selected organism.",
          type = "warning"
        )
        return()
      }
      shiny::updateSelectInput(session, "datatype", selected = "RNA-seq")
      shiny::updateTextAreaInput(session, "features",
        value = paste(features, collapse = "\n")
      )
    })

    result <- shiny::eventReactive(input$convert, {
      probes <- parse_feature_list(input$features)
      shiny::validate(shiny::need(
        length(probes) > 0,
        "Please paste at least one gene/feature ID."
      ))

      annot <- shiny::withProgress(
        message = "Converting gene/feature IDs...",
        playbase::getProbeAnnotation(
          organism = input$organism,
          probes = probes,
          datatype = input$datatype
        )
      )

      shiny::validate(shiny::need(
        !is.null(annot),
        "Could not annotate these IDs for the selected organism."
      ))

      annot
    })

    output$table_area <- shiny::renderUI({
      df <- tryCatch(result(), error = function(e) NULL)
      if (is.null(df)) {
        shiny::div(
          style = paste(
            "display: flex; flex-direction: column; align-items: center;",
            "justify-content: center; height: 100%; color: #aaa;",
            "font-size: 14px; text-align: center;"
          ),
          shiny::icon("table", class = "fa-2x", style = "margin-bottom: 10px;"),
          "Paste gene/feature IDs and click Convert to see results here."
        )
      } else {
        DT::dataTableOutput(session$ns("table"), height = "100%")
      }
    })

    output$table <- DT::renderDataTable({
      df <- result()
      df$gene_name <- NULL
      
      DT::datatable(
        data = df,
#        class = "compact hover",
        rownames = FALSE,
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = 'none',
        options = list(
          dom = "lfrtip",
          scroller = TRUE,
          scrollX = TRUE,
          scrollY = "calc(100vh - 60px)",
          scrollResize = TRUE
#         deferRender = TRUE,
#         autoWidth = TRUE
        )
      ) %>% DT::formatStyle(
        columns = 0,
        target = "row",
        ## fontSize = "14px",
        lineHeight = "90%"
      ) 

    })

    output$download_ui <- shiny::renderUI({
      df <- tryCatch(result(), error = function(e) NULL)
      shiny::downloadButton(session$ns("download"), "Download CSV",
        enabled = !is.null(df), style = "width: 100%;"
      )
    })

    output$download <- shiny::downloadHandler(
      filename = function() "feature_annotation.csv",
      content = function(file) {
        utils::write.csv(result(), file, row.names = FALSE)
      }
    )
  })
}
