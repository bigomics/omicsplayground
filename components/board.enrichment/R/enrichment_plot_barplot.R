##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

enrichment_plot_barplot_ui <- function(id, height, width) {
  ns <- shiny::NS(id)

  info_text <- "An enrichment barplot per sample group for the gene set that is selected from the enrichment analysis Table <code>I</code>. Samples can be ungrouped in the barplot by selecting <code>ungroup samples</code> from the plot <i>Settings</i>."

  options <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("gs_ungroup1"), "ungroup samples", FALSE
      ),
      "Ungroup samples in the plot",
      placement = "top",
      options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = "Enrichment barplot",
    label = "b",
    info.text = info_text,
    options = options,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

enrichment_plot_barplot_server <- function(id,
                                           inputData,
                                           gset_selected,
                                           gs_contrast,
                                           subplot.MAR,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    subplot_barplot.RENDER <- shiny::reactive({
      par(mfrow = c(1, 1), mgp = c(1.8, 0.8, 0), oma = c(0, 0, 0, 0.4))
      par(mar = subplot.MAR)

      ngs <- inputData()
      shiny::req(ngs)


      gset <- rownames(ngs$gsetX)[1]
      gset <- gset_selected()
      if (is.null(gset) || length(gset) == 0) {
        return(NULL)
      }
      gset <- gset[1]
      if (!gset %in% rownames(ngs$gsetX)) {
        return(NULL)
      }

      comp0 <- colnames(ngs$model.parameters$contr.matrix)[1]
      comp0 <- gs_contrast()

      grouped <- TRUE
      grouped <- FALSE
      grouped <- !input$gs_ungroup1
      has.design <- !is.null(ngs$model.parameters$design)
      collapse.others <- ifelse(has.design, FALSE, TRUE)
      ## collapse.others=TRUE

      ngrp <- length(unique(ngs$samples$group))
      srt <- ifelse(!grouped || ngrp > 4, 30, 0)
      if (!grouped && ncol(ngs$X) > 15) srt <- 60
      pgx.plotExpression(
        ngs, gset,
        comp = comp0, logscale = TRUE, level = "geneset",
        collapse.others = collapse.others, grouped = grouped,
        cex = 1.1, srt = srt, main = "", ylab = "enrichment (avg logFC)"
      )
      title(breakstring(gset, 42, 80), cex.main = 0.85)
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = subplot_barplot.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 100),
      add.watermark = watermark
    )
  })
}
