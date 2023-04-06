##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_top_enrich_gsets_ui <- function(id, height, width) {
  ns <- shiny::NS(id)

  info_text <- "This plot shows the <strong>top enriched</strong> gene sets for the selected comparison in the <code>Contrast</code> settings. Black vertical bars indicate the rank of genes in the gene set in the sorted list metric. The green curve corresponds to the 'running statistics' of the enrichment score (ES). The more the green ES curve is shifted to the upper left of the graph, the more the gene set is enriched in the first group. Conversely, a shift of the ES curve to the lower right, corresponds to more enrichment in the second group."

  PlotModuleUI(
    ns("plot"),
    title = "Top enriched gene sets",
    label = "a",
    info.text = info_text,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

enrichment_plot_top_enrich_gsets_server <- function(id,
                                                    pgx,
                                                    getFilteredGeneSetTable,
                                                    gs_contrast,
                                                    gseatable,
                                                    watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plotTopEnriched <- function(pgx, rpt, comp, ntop, rowcol) {
      if (is.null(pgx)) {
        return(NULL)
      }

      gx.meta <- pgx$gx.meta$meta[[comp]]
      rnk0 <- gx.meta$meta.fx
      names(rnk0) <- pgx$genes[rownames(gx.meta), "gene_name"]
      rnk0 <- rnk0 - mean(rnk0, na.rm = TRUE) ## scaling/centering should be done in calculation...
      fx.col <- grep("score|fx|fc|sign|NES|logFC", colnames(rpt))[1]
      qv.col <- grep("meta.q|q$", colnames(rpt))[1]
      fx <- rpt[, fx.col]
      qv <- rpt[, qv.col]
      names(qv) <- names(fx) <- rownames(rpt)
      top <- rownames(rpt)
      top <- setdiff(top, c(NA, "NA"))
      if (is.null(top) || is.na(top[1])) {
        return(NULL)
      }

      par(mfrow = rowcol)
      if (ntop == 1) {
        par(mar = c(1, 6, 2, 6), mgp = c(1.6, 0.6, 0), oma = c(0.1, 1, 0, 0.1))
      } else {
        par(mar = c(0.2, 1.8, 2.3, 0.1), mgp = c(1.6, 0.6, 0), oma = c(0.1, 1, 0, 0.1))
      }

      for (i in 1:ntop) {
        gs <- top[i]
        if (i > length(top) || is.na(gs)) {
          frame()
        } else {
          genes <- names(which(pgx$GMT[, gs] != 0))
          genes <- toupper(genes)
          names(rnk0) <- toupper(names(rnk0))
          ylab <- ""
          if (i %% rowcol[2] == 1) ylab <- "Rank metric"
          xlab <- ""
          gs1 <- playbase::breakstring(gs, 28, 50, force = FALSE)
          if (ntop == 1) {
            gs1 <- playbase::breakstring(gs, 100, 200, force = FALSE)
            xlab <- "Rank in ordered dataset"
            ylab <- "Rank metric"
          }
          playbase::gsea.enplot(rnk0, genes,
            names = NULL, ## main=gs,
            main = gs1, xlab = xlab, ylab = ylab,
            lab.line = c(0, 1.8), cex.lab = 0.75,
            cex.main = 0.78, len.main = 200
          )
          qv1 <- formatC(qv[gs], format = "e", digits = 2)
          legend("topright", paste("q=", qv1), bty = "n", cex = 0.85)
        }
      }
    }

    topEnriched.RENDER <- shiny::reactive({
      shiny::req(pgx)
      rpt <- getFilteredGeneSetTable()
      shiny::req(rpt, gs_contrast())
      if (is.null(rpt)) {
        return(NULL)
      }

      comp <- 1
      comp <- gs_contrast()
      if (!(comp %in% names(pgx$gx.meta$meta))) {
        return(NULL)
      }

      ## selected
      sel <- as.integer(gseatable$rows_selected())
      sel.gs <- NULL
      if (!is.null(sel) && length(sel) > 0) sel.gs <- rownames(rpt)[sel]

      ii <- gseatable$rows_selected()
      jj <- gseatable$rows_current()
      shiny::req(jj)

      if (nrow(rpt) == 0) {
        return(NULL)
      }

      ## ENPLOT TYPE
      if (length(ii) > 0) {
        itop <- ii[1]
      } else {
        itop <- head(jj, 15)
      }
      if (length(itop) == 1) {
        plotTopEnriched(pgx, rpt[itop, , drop = FALSE], comp = comp, ntop = 1, rowcol = c(1, 1))
      } else {
        plotTopEnriched(pgx, rpt[itop, , drop = FALSE], comp = comp, ntop = 15, rowcol = c(3, 5))
      }
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = topEnriched.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(90, 120),
      add.watermark = watermark
    )
  })
}
