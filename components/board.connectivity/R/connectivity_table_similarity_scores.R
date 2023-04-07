##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

connectivity_table_similarity_scores_ui <- function(
  id,
  title, 
  info.text,
  caption,
  width,
  height,
  label="") {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    width = width,
    caption = caption,
    height = height,
    title = title,
    label = label
  )
}

connectivity_table_similarity_scores_server <- function(id,
                                                        getConnectivityScores,
                                                        columns,
                                                        height) {
  moduleServer(id, function(input, output, session) {
    connectivityScoreTable.RENDER <- shiny::reactive({
      df <- getConnectivityScores()
      shiny::req(df)

      ##kk <- c("pathway", "score", "rho", "NES", "padj", "leadingEdge")
      kk <- intersect(columns, colnames(df))
      df <- df[, kk]
      df <- df[abs(df$score) > 0, , drop = FALSE]

      df$pathway <- playbase::shortstring(df$pathway, 100)

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
          scrollY = height,
          scroller = TRUE, deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numcols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("score",
          background = playbase::color_from_middle(
            df[, "score"], "lightblue", "#f5aeae"
          ),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    connectivityScoreTable.RENDER_modal <- shiny::reactive({
      dt <- connectivityScoreTable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    connectivityScoreTable <- TableModuleServer(
      "datasets",
      func = connectivityScoreTable.RENDER,
      func2 = connectivityScoreTable.RENDER_modal,
      selector = "single"
    )

    return(connectivityScoreTable)
  })
}
