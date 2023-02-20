##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

wgcna_table_enrichment_ui <- function(id, width, height) {
  ns <- shiny::NS(id)

  info_text <- "In this table, users can check mean expression values of features across the conditions for the selected genes."

  TableModuleUI(
    ns("datasets"),
    info.text = info_text,
    width = width,
    height = height,
    title = "Module enrichment",
    label = "e"
  )

}

wgcna_table_enrichment_server <- function(id,
                                          enrich_table) {
  moduleServer(id, function(input, output, session) {
    enrichTable.RENDER <- shiny::reactive({
      df <- enrich_table()
      numeric.cols <- grep("score|value|ratio", colnames(df))

      DT::datatable(
        df,
        rownames = FALSE,
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE,
          scrollY = "70vh",
          scroller = TRUE, deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    })

    enrichTable.RENDER_modal <- shiny::reactive({
      dt <- enrichTable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    enrichTable_module <- TableModuleServer(
      "datasets",
      func = enrichTable.RENDER,
      func2 = enrichTable.RENDER_modal,
      selector = "none"
    )

    return(enrichTable_module)
  })
}
