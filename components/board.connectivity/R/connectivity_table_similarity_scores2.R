##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

connectivity_table_similarity_scores2_ui <- function(id, width, height) {
  ns <- shiny::NS(id)

  connectivityScoreTable_info <- "<b>Similarity scores.</b> Normalized enrichment scores (NES) and Pearson correlation (rho) of reference profiles with respect to the currently selected contrast. The top 100 up/down genes are considered for the calculation of rho or NES. The score is calculated as rho^2*NES. "

  TableModuleUI(
    ns("datasets"),
    info.text = connectivityScoreTable_info,
    width = width,
    height = height,
    title = "Similarity scores",
    label = "b"
  )
}

connectivity_table_similarity_scores2_server <- function(id,
                                                         getConnectivityScores) {
  moduleServer(id, function(input, output, session) {
    connectivityScoreTable2.RENDER <- shiny::reactive({
      df <- getConnectivityScores()
      if (is.null(df)) {
        return(NULL)
      }

      kk <- c("pathway", "score", "rho", "NES", "padj", "size", "leadingEdge")
      kk <- c("score", "pathway", "rho", "NES", "padj")
      kk <- intersect(kk, colnames(df))
      df <- df[, kk]
      df <- df[abs(df$score) > 0, , drop = FALSE]

      df$pathway <- shortstring(df$pathway, 110)
      df$leadingEdge <- NULL

      colnames(df) <- sub("pathway", "dataset/contrast", colnames(df))
      score.col <- which(colnames(df) == "score")
      numcols <- c("score", "pval", "padj", "NES.q", "ES", "NES", "rho", "R2")
      numcols <- intersect(numcols, colnames(df))

      DT::datatable(df,
        rownames = FALSE,
        class = "compact cell-border stripe hover",
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          pageLength = 99999,
          scrollX = TRUE,
          scrollY = "55vh",
          scroller = TRUE, deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numcols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("score",
          background = color_from_middle(
            df[, "score"], "lightblue", "#f5aeae"
          ),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    connectivityScoreTable2.RENDER_modal <- shiny::reactive({
      dt <- connectivityScoreTable2.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    connectivityScoreTable <- TableModuleServer(
      "datasets",
      func = connectivityScoreTable2.RENDER,
      func2 = connectivityScoreTable2.RENDER_modal,
      selector = "none"
    )

    return(connectivityScoreTable)
  })
}
