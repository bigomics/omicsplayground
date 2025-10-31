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
  width
) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::selectInput( ns("plottype"), "Plot type:",
      choices = c("geneset score","gene frequency"),
      selected = "gene frequency"
    )

  )

  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

wgcna_plot_enrichment_server <- function(id,
                                         enrichTable,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    enrichPlot.RENDER <- function() {
      df <- enrichTable$data()
      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }
      ii <- enrichTable$rows_all()
      shiny::req(ii)
      df <- df[ii, , drop = FALSE]
      df <- head(df, 20)

      if(input$plottype == "gene frequency") {
        genes <- unlist(strsplit(df$genes,split="\\|"))
        gtable <- sort(table(genes),decreasing=TRUE)
        par(mar = c(8,4,2,0))
        barplot(
          head(gtable,35),
          las = 3,
          cex.names = 0.9,
          xlab = "",
          ylab = "Frequency in top"
        )
      }

      if(input$plottype == "geneset score") {
        gs.top <- df$geneset
        xlim0 <- c(0, max(df$score))
        col1 <- c("grey90", "#f5bfbf")[1 + 1 * (df$q.value < 0.05)]
        par(mar = c(4.5, 1, 1, 1))
        barplot(rev(df$score),
          horiz = TRUE, width = 0.8, space = 0.25, xlim = xlim0,
          border = NA, col = rev(col1), xlab = "enrichment score"
        )
        text(0, (nrow(df):1) - 0.48, gs.top, adj = 0, pos = 4, cex = 0.8)
      }
      
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
