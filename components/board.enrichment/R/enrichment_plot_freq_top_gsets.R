##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

enrichment_plot_freq_top_gsets_ui <- function(id, height, width) {
  ns <- shiny::NS(id)

  info_text <- "<strong>Gene frequency.</strong> The plot shows the number of times a gene is present in the top-N genesets sorted by frequency. Genes that are frequently shared among the top enriched gene sets may suggest driver genes."

  topEnrichedFreq.opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(ns("gs_enrichfreq_ntop"), "Number of top sets",
        c(5, 10, 15),
        inline = TRUE, selected = 15
      ),
      "Number of top genesets to consider for counting the gene frequency."
    ),
    withTooltip(
      shiny::checkboxInput(
        ns("gs_enrichfreq_gsetweight"),
        "Weight by geneset size", TRUE
      ),
      "Weight by (inverse) gene set size."
    ),
    withTooltip(
      shiny::checkboxInput(
        ns("gs_enrichfreq_fcweight"),
        "Weight by FC", TRUE
      ),
      "Weight by fold-change of current contrast."
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = "Frequency in top gene sets",
    label = "b",
    info.text = info_text,
    options = topEnrichedFreq.opts,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "csv")
  )
}

enrichment_plot_freq_top_gsets_server <- function(id,
                                                  inputData,
                                                  getFilteredGeneSetTable,
                                                  gs_contrast,
                                                  gseatable,
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plotEnrichFreq <- function(ngs, rpt, ntop, ngenes, gset.weight, fcweight) {
      fx.col <- grep("score|fx|fc|sign|NES|logFC", colnames(rpt))[1]
      fx <- rpt[, fx.col]
      names(fx) <- rownames(rpt)

      top <- rownames(rpt)
      top <- head(top, ntop)
      if (!all(top %in% colnames(ngs$GMT))) {
        return(NULL)
      }

      F <- 1 * (ngs$GMT[, top, drop = FALSE] > 0)
      F <- as.matrix(F)
      wt <- FALSE
      if (gset.weight) {
        F <- Matrix::t(Matrix::t(F) / Matrix::colSums(F, na.rm = TRUE))
        wt <- TRUE
      }
      F <- Matrix::t(Matrix::t(F) * sign(fx[top]))
      if (fcweight) {
        F <- Matrix::t(Matrix::t(F) * abs(fx[top]))
        wt <- TRUE
      }
      F <- head(F[order(-Matrix::rowSums(abs(F))), , drop = FALSE], ngenes)
      F <- F[order(-Matrix::rowSums(F)), , drop = FALSE]

      sel.zero <- which(Matrix::rowSums(abs(F)) < 1e-4)
      if (length(sel.zero)) rownames(F)[sel.zero] <- ""

      par(mfrow = c(1, 1), mar = c(6, 4, 2, 0.5), mgp = c(2.2, 0.8, 0))
      col1 <- grey.colors(ncol(F), start = 0.15)
      ylab <- ifelse(wt, "weighted frequency", "frequency")
      barplot(t(F),
        beside = FALSE, las = 3, cex.names = 0.90, col = col1,
        ylab = ylab
      )
    }

    plot_data <- shiny::reactive({
      ngs <- inputData()

      rpt <- getFilteredGeneSetTable()
      shiny::req(ngs, rpt, gs_contrast())

      comp <- gs_contrast()
      if (is.null(comp)) {
        return(NULL)
      }
      if (!(comp %in% names(ngs$gx.meta$meta))) {
        return(NULL)
      }

      ## filter on active rows (using search)
      ii <- gseatable$rows_current()
      rpt <- rpt[ii, , drop = FALSE]
      if (nrow(rpt) == 0) {
        return(NULL)
      }
      ntop <- as.integer(input$gs_enrichfreq_ntop)
      gset.weight <- input$gs_enrichfreq_gsetweight
      fcweight <- input$gs_enrichfreq_fcweight
      return(
        list(
          ngs,
          rpt,
          ntop,
          gset.weight,
          fcweight
        )
      )
    })

    topEnrichedFreq.RENDER <- shiny::reactive({
      dt <- plot_data()
      ngs <- dt[[1]]
      rpt <- dt[[2]]
      ntop <- dt[[3]]
      gset.weight <- dt[[4]]
      fcweight <- dt[[5]]
      plotEnrichFreq(ngs, rpt, ntop = ntop, ngenes = 35, gset.weight, fcweight)
      p <- grDevices::recordPlot()
      p
    })

    topEnrichedFreq.RENDER2 <- shiny::reactive({
      dt <- plot_data()
      ngs <- dt[[1]]
      rpt <- dt[[2]]
      ntop <- dt[[3]]
      gset.weight <- dt[[4]]
      fcweight <- dt[[5]]
      plotEnrichFreq(ngs, rpt, ntop = ntop, ngenes = 60, gset.weight, fcweight)
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = topEnrichedFreq.RENDER,
      func2 = topEnrichedFreq.RENDER2,
      pdf.width = 5, pdf.height = 5,
      res = c(68, 100),
      csvFunc = plot_data,
      add.watermark = watermark
    )
  })
}
