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
      shiny::checkboxInput(ns("FoldchangeHeatmap_allfc"), "show all contrasts", TRUE),
      "Show all contrasts or just the selected ones."
    ),
    withTooltip(
      shiny::checkboxInput(ns("FoldchangeHeatmap_cluster"), "cluster genes", FALSE),
      "Cluster genes (columns)."
    ),
    withTooltip(
      shiny::radioButtons(ns("FoldchangeHeatmap_annotype"), "annotation type",
        c("boxplot", "barplot"),
        inline = TRUE
      ),
      "Plot type of column annotation."
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
                                      watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      if (input$FoldchangeHeatmap_allfc) {
        F <- getFoldChangeMatrix()$fc
      } else {
        F <- getActiveFoldChangeMatrix()$fc
      }

      F <- F[order(-rowMeans(F**2)), ]
      F <- F[order(-abs(rowMeans(F))), ]

      F1 <- head(F, 80)
      F1 <- F1[order(rowMeans(F1)), ]
      F1
    })

    FoldchangeHeatmap.PLOT <- function() {
      F1 <- plot_data()
      bh <- 5
      mh <- 6
      bm <- 4
      cclust <- TRUE
      mh <- min(max(ncol(F1) * 0.35, 0.8), 8)
      cclust <- input$FoldchangeHeatmap_cluster
      if (level == "geneset") {
        bh <- 3
      }
      bm <- 5 - mh ## bottom margin
      at <- input$FoldchangeHeatmap_annotype

      par(mfrow = c(1, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 3, 0))

      plt <- grid::grid.grabExpr({
        frame()
        playbase::heatmapWithAnnot(
          F1,
          anno.type = at,
          bar.height = bh, map.height = mh,
          mar = c(bm, 0.5, 0.5, 1),
          cluster_columns = cclust,
          inset = c(0.01, 0.01)
        )
      })

      # plt
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
