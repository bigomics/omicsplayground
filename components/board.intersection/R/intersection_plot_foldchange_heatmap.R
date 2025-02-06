##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

foldchange_heatmap_ui <- function(
    id,
    title,
    info.text,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  FoldchangeHeatmap.opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(ns("allfc"), "show all contrasts", FALSE),
      "Show all contrasts or just the selected ones."
    ),
    withTooltip(
      shiny::checkboxInput(ns("cluster"), tspan("cluster genes"), FALSE),
      "Cluster genes (columns)."
    )
  )

  PlotModuleUI(
    ns("FoldchangeHeatmap"),
    title = title,
    label = "a",
    plotlib = "grid",
    info.text = info.text,
    caption = caption,
    options = FoldchangeHeatmap.opts,
    download.fmt = c("png", "pdf", "csv"),
    height = height,
    width = width
  )
}


foldchange_heatmap_server <- function(id,
                                      getFoldChangeMatrix,
                                      getActiveFoldChangeMatrix,
                                      pgx,
                                      level,
                                      input_comparisons,
                                      watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      if (input$allfc) {
        FC <- getFoldChangeMatrix()$fc
      } else {
        FC <- getActiveFoldChangeMatrix()$fc
        ## FC <- FC[, input_comparisons(), drop = FALSE]
      }

      # check that dim(F)[2] >0, in case user has not selected any comparisons
      validate(
        need(
          dim(FC)[2] > 0,
          "No comparisons selected. Please select a comparison on the settings sidebar."
        )
      )

      FC <- FC[order(-rowMeans(FC**2, na.rm = TRUE)), , drop = FALSE]
      FC <- FC[order(-abs(rowMeans(FC, na.rm = TRUE))), , drop = FALSE]

      F1 <- head(FC, 80)
      F1 <- F1[order(rowMeans(F1, na.rm = TRUE)), , drop = FALSE]
      F1
    })

    FoldchangeHeatmap.PLOT <- function() {
      F1 <- plot_data()
      shiny::req(F1)
      bh <- 5
      mh <- 6
      bm <- 4
      cclust <- TRUE
      mh <- min(max(ncol(F1) * 0.35, 0.8), 8)
      cclust <- input$cluster
      if (level == "geneset") {
        bh <- 3
      }
      bm <- 5 - mh ## bottom margin
      par(mfrow = c(1, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 3, 0))

      plt <- grid::grid.grabExpr({
        frame()
        playbase::heatmapWithAnnot(
          F1,
          anno.type = "boxplot",
          bar.height = bh,
          map.height = mh,
          mar = c(bm, 0.5, 0.5, 1),
          cluster_columns = cclust,
          inset = c(0.01, 0.01)
        )
      })

      grid::grid.draw(plt)
    }

    PlotModuleServer(
      "FoldchangeHeatmap",
      func = FoldchangeHeatmap.PLOT,
      csvFunc = plot_data,
      plotlib = "base",
      res = c(90, 110),
      pdf.width = 5, pdf.height = 5,
      add.watermark = watermark
    )
  })
}
