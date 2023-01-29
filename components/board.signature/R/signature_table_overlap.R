##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

signature_table_overlap_ui <- function(id) {
  ns <- shiny::NS(id)

  tableWidget(ns("table"))
}

signature_table_overlap_server <- function(id,
                                           getOverlapTable,
                                           fullH,
                                           tabH) {
  moduleServer(id, function(input, output, session) {
    overlapTable.RENDER <- shiny::reactive({
      df <- getOverlapTable()
      shiny::req(df)

      df$geneset <- wrapHyperLink(df$geneset, df$geneset)

      numeric.cols <- which(sapply(df, is.numeric))
      numeric.cols <- intersect(c("p.fisher", "q.fisher"), colnames(df))

      DT::datatable(df,
        class = "compact cell-border stripe",
        rownames = FALSE, escape = c(-1, -2),
        extensions = c("Scroller"),
        selection = "none",
        fillContainer = TRUE,
        options = list(
          dom = "frtip",
          scrollX = TRUE, scrollY = tabH, scroller = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("score",
          background = color_from_middle(df$score, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    info.text <- "Under the <strong>Overlap/similarity tab</strong>, users can find the similarity of their gene list with all the gene sets and pathways in the platform, including statistics such as the total number of genes in the gene set (K), the number of intersecting genes between the list and the gene set (k), the overlapping ratio of k/K, logarithm of the  odds ratio (log.OR), as well as the p and q values by the Fisher’s test for the overlap test."

    overlapTable <- shiny::callModule(
      tableModule,
      id = "table",
      func = overlapTable.RENDER,
      title = tags$div(
        HTML('<span class="module-label">(b)</span>Overlap with other signatures')
      ),
      info.text = info.text,
      height = 0.4 * fullH
    )
    return(overlapTable)
  })
}
