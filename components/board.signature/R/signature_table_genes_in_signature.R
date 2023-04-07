##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

signature_table_genes_in_signature_ui <- function(
  id,
  title,
  info.text,
  caption,
  width,
  height) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    width = width,
    height = height,
    title = title,
    caption, caption,
    label = "b"
  )
}

signature_table_genes_in_signature_server <- function(id,
                                                      getEnrichmentGeneTable
                                                      ) {
  moduleServer(id, function(input, output, session) {

    enrichmentGeneTable.RENDER <- shiny::reactive({
      df <- getEnrichmentGeneTable()
      if (is.null(df)) {
        shiny::validate(shiny::need(!is.null(df), "Please select a contrast"))
        return(NULL)
      }

      color_fx <- as.numeric(df[, 3:ncol(df)])
      color_fx[is.na(color_fx)] <- 0 ## yikes...

      numeric.cols <- colnames(df)[3:ncol(df)]

      DT::datatable(df,
        class = "compact cell-border stripe",
        rownames = FALSE,
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        fillContainer = TRUE,
        options = list(
          dom = "lrftip",
          scrollX = TRUE, scrollY = "30vh", scroller = TRUE,
          deferRender = FALSE
        )
      ) %>% ## end of options.list
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(
          numeric.cols,
          background = playbase::color_from_middle(color_fx, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    enrichmentGeneTable.RENDER_modal <- shiny::reactive({
      dt <- enrichmentGeneTable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    enrichmentGeneTable <- TableModuleServer(
      "datasets",
      func = enrichmentGeneTable.RENDER,
      func2 = enrichmentGeneTable.RENDER_modal,
      selector = "single"
    )
    return(enrichmentGeneTable)
  })
}
