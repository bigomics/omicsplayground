##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

wordcloud_table_enrichment_ui <- function(id) {
  ns <- shiny::NS(id)

  tableWidget(ns("wordcloud_enrichmentTable"))
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
        selection = list(mode = "single", target = "row", selected = 1),
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE, scrollY = 200,
          scroller = TRUE, deferRender = TRUE
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

    wordcloud_enrichmentTable_info <-
      "<b>Keyword enrichment table.</b> This table shows the keyword enrichment statistics for the selected contrast. The enrichment is calculated using GSEA for occurance of the keywork in the ordered list of gene set descriptions."

    wordcloud_enrichmentTable <- shiny::callModule(
      tableModule,
      id = "wordcloud_enrichmentTable",
      func = wordcloud_enrichmentTable.RENDER,
      info.text = wordcloud_enrichmentTable_info,
      title = tags$div(
        HTML('<span class="module-label">(e)</span>Enrichment table')
      ),
      height = c(270, 700)
    )

    return(wordcloud_enrichmentTable)
  })
}
