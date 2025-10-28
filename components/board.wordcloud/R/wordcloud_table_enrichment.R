##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wordcloud_table_enrichment_ui <- function(
  id,
  title,
  info.text,
  caption,
  width,
  height
) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    width = width,
    height = height,
    title = title,
    caption = caption,
    label = "d"
  )
}

wordcloud_table_enrichment_server <- function(id,
                                              getCurrentWordEnrichment) {
  moduleServer(id, function(input, output, session) {
    wordcloud_enrichmentTable.RENDER <- shiny::reactive({
      df <- getCurrentWordEnrichment()
      shiny::req(df)
      df <- df[, c("word", "pval", "padj", "ES", "NES", "size")]

      numeric.cols <- colnames(df)[which(sapply(df, is.numeric))]
      numeric.cols
      tbl <- DT::datatable(
        df,
        rownames = FALSE,
        class = "compact cell-border stripe hover",
        extensions = c("Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = 1),
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE,
          scrollY = "25vh",
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("NES",
          background = color_from_middle(df[, "NES"], "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
      return(tbl)
    })

    wordcloud_enrichmentTable.RENDER_modal <- shiny::reactive({
      dt <- wordcloud_enrichmentTable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    wordcloud_enrichmentTable <- TableModuleServer(
      "datasets",
      func = wordcloud_enrichmentTable.RENDER,
      func2 = wordcloud_enrichmentTable.RENDER_modal,
      selector = "single"
    )

    return(wordcloud_enrichmentTable)
  })
}
