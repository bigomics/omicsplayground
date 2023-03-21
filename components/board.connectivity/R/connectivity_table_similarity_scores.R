##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

connectivity_table_similarity_scores_ui <- function(id, width, height) {
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

connectivity_table_similarity_scores_server <- function(id,
                                                        getConnectivityScores,
                                                        cmap_sigdb) {
  moduleServer(id, function(input, output, session) {
    connectivityScoreTable.RENDER <- shiny::reactive({
      df <- getConnectivityScores()
      shiny::req(df)

      kk <- c("pathway", "score", "rho", "NES", "padj", "leadingEdge")
      kk <- intersect(kk, colnames(df))
      df <- df[, kk]
      df <- df[abs(df$score) > 0, , drop = FALSE]

      ## --------- temporarily add LINCS descriptive name !!!!!!!!!!!!!! -----------------
      if (DEV && cmap_sigdb() == "sigdb-lincs.h5" && !is.null(PERTINFO)) {
        dd <- sub("\\|.*", "", df$pathway)
        pert_iname <- PERTINFO[match(dd, rownames(PERTINFO)), "pert_iname"]
        df$pathway <- paste0(df$pathway, " (", pert_iname, ")")
      }
      ## ---------- temporarily add LINCS descriptive name !!!!!!!!!!!!!! -----------------

      ## colnames(df) <- sub("padj","NES.q",colnames(df))
      df$leadingEdge <- shortstring(sapply(df$leadingEdge, paste, collapse = ","), 40)
      df$pathway <- shortstring(df$pathway, 100)
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
          scrollY = "25vh",
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
