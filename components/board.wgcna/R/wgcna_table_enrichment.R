##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_table_enrichment_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    width,
    height) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    title = title,
    info.text = info.text,
    caption = caption,
    width = width,
    height = height,
    label = label
  )
}

wgcna_table_enrichment_server <- function(id,
                                          enrich_table) {
  moduleServer(id, function(input, output, session) {

    RENDER <- shiny::reactive({
      df <- enrich_table()
      numeric.cols <- grep("score|value|ratio", colnames(df))

      DT::datatable(
        df,
        rownames = FALSE,
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE,
          scrollY = "70vh",
          scrollResize = TRUE,
          scroller = TRUE, deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "10px", lineHeight = "70%")
    })

    RENDER_modal <- shiny::reactive({
      dt <- RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    tablemodule <- TableModuleServer(
      "datasets",
      func = RENDER,
      func2 = RENDER_modal,
      selector = "none"
    )

    return(tablemodule)
  })
}
