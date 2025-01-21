##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_enrichment_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height,
    width) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

wgcna_plot_enrichment_server <- function(id,
                                         enrich_table,
                                         enrichTable_module,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    enrichPlot.RENDER <- function() {
      df <- enrich_table()
      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }

      ii <- enrichTable_module$rows_all()
      shiny::req(ii)
      df <- df[ii, , drop = FALSE]
      df <- head(df, 20)      
      gs.top <- df$geneset
      xlim0 <- c(0, max(df$score))
      col1 <- c("lightskyblue1", "lightpink")[1 + 1 * (df$q.value < 0.05)]
      par(mar = c(4.5, 1, 1, 1))
      barplot(rev(df$score),
        horiz = TRUE, width = 0.8, space = 0.25, xlim = xlim0,
        border = NA, col = rev(col1), xlab = "score  (odd.ratio * -log10p)"
      )
      text(0, (nrow(df):1) - 0.48, gs.top, adj = 0, pos = 4, cex = 0.8)
    }

    PlotModuleServer(
      "plot",
      func = enrichPlot.RENDER,
      pdf.width = 10, pdf.height = 5,
      res = c(80, 110),
      add.watermark = watermark
    )
  })
}
