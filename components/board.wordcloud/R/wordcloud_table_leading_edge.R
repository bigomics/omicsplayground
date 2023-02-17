##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

wordcloud_table_leading_edge_ui <- function(id, width, height) {
  ns <- shiny::NS(id)

  info_text <- "Keyword leading edge table."

  TableModuleUI(
    ns("datasets"),
    info.text = info_text,
    width = width,
    height = height,
    title = "Leading-edge table",
    label = "e"
  )

}

wordcloud_table_leading_edge_server <- function(id,
                                                pgx,
                                                wc_contrast,
                                                wordcloud_enrichmentTable,
                                                getCurrentWordEnrichment) {
  moduleServer(id, function(input, output, session) {
    wordcloud_leadingEdgeTable.RENDER <- shiny::reactive({
      shiny::req(pgx$gset.meta, wc_contrast())

      df <- getCurrentWordEnrichment()
      sel.row <- wordcloud_enrichmentTable$rows_selected()
      shiny::req(df, sel.row)

      ee <- unlist(df$leadingEdge[sel.row])
      ee <- strsplit(ee, split = "//")[[1]]

      fx <- pgx$gset.meta$meta[[wc_contrast()]][ee, "meta.fx"]
      names(fx) <- ee
      df <- data.frame("leading.edge" = ee, fx = fx)
      df <- df[order(-abs(df$fx)), ]
      rownames(df) <- ee

      numeric.cols <- colnames(df)[which(sapply(df, is.numeric))]

      df$leading.edge <- wrapHyperLink(df$leading.edge, df$leading.edge) ## add link

      tbl <- DT::datatable(df,
        rownames = FALSE, escape = c(-1, -2),
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
        DT::formatStyle("fx",
          background = color_from_middle(df[, "fx"], "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
      return(tbl)
    })

    wordcloud_leadingEdgeTable <- TableModuleServer(
      "datasets",
      func = wordcloud_leadingEdgeTable.RENDER,
      selector = "none"
    )

    return(wordcloud_leadingEdgeTable)
  })
}
