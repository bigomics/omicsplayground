##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

enrichment_plot_volcano_ui <- function(id, height, width) {
  ns <- shiny::NS(id)

  info_text <- "<b>Volcano plot.</b> Volcano-plot showing significance versus fold-change on the y and x axes, respectively. Genes in the gene set that is selected from the enrichment analysis <b>Table I</b> are highlighted in blue."

  PlotModuleUI(
    ns("plot"),
    title = "Volcano plot",
    label = "a",
    info.text = info_text,
    height = height,
    width = width,
    plotlib = "base",
    plotlib2 = "plotly",
    download.fmt = c("png", "pdf")
  )
}

enrichment_plot_volcano_server <- function(id,
                                           inputData,
                                           gs_contrast,
                                           selected_gxmethods,
                                           gset_selected,
                                           gs_fdr,
                                           gs_lfc,
                                           subplot.MAR,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    subplot_volcano.RENDER <- shiny::reactive({
      par(mfrow = c(1, 1), mgp = c(1.2, 0.4, 0), oma = c(0, 0, 0, 0.4))
      par(mar = subplot.MAR)

      ngs <- inputData()
      shiny::req(ngs)

      comp <- 1
      gs <- 1
      comp <- gs_contrast()
      ngs <- inputData()
      shiny::req(ngs)

      gxmethods <- selected_gxmethods() ## from module-expression
      shiny::req(gxmethods)

      gx.meta <- ngs$gx.meta$meta[[comp]]
      meta.q <- apply(gx.meta$q[, gxmethods, drop = FALSE], 1, max) ## max q-value
      limma1 <- data.frame(meta.fx = gx.meta$meta.fx, meta.q = meta.q)
      gx.annot <- ngs$genes[rownames(gx.meta), c("gene_name", "gene_title")]
      limma <- cbind(gx.annot, limma1)

      gs <- gset_selected()
      if (is.null(gs) || length(gs) == 0) {
        frame()
        text(0.5, 0.5, "Please select a geneset", col = "grey50")
        return()
      }
      gs <- gs[1]

      gset <- getGSETS(gs)[[1]]
      dbg("[subplot_volcano.RENDER] head.gset = ", head(gset, 5))
      jj <- match(toupper(gset), toupper(limma$gene_name))
      sel.genes <- setdiff(limma$gene_name[jj], c(NA, "", " "))

      dbg("[subplot_volcano.RENDER] head.sel.genes = ", head(sel.genes, 5))

      fdr <- 1
      fdr <- as.numeric(gs_fdr())
      fc.genes <- as.character(limma[, grep("^gene$|gene_name", colnames(limma))])
      fx <- limma[, grep("logFC|meta.fx|fc", colnames(limma))[1]]
      qval <- limma[, grep("^q|adj.P.Val|meta.q|qval|padj", colnames(limma))[1]]

      qval <- pmax(qval, 1e-12) ## prevent q=0
      qval[which(is.na(qval))] <- 1
      xlim <- c(-1, 1) * max(abs(fx), na.rm = TRUE)
      ylim <- c(0, 12)
      ylim <- c(0, max(12, 1.1 * max(-log10(qval), na.rm = TRUE)))

      lfc <- 0.20
      lfc <- as.numeric(gs_lfc())

      gx.volcanoPlot.XY(
        x = fx, pv = qval, gene = fc.genes,
        render = "canvas", n = 5000, nlab = 10,
        xlim = xlim, ylim = ylim, ## hi.col="#222222",
        use.fdr = TRUE, p.sig = fdr, lfc = lfc,
        cex = 0.9, lab.cex = 1.3,
        cex.main = 0.8, cex.axis = 0.9,
        xlab = "fold change (log2)",
        ylab = "significance (log10q)",
        highlight = sel.genes
      )
      gs <- breakstring(gs, 50)
      title(gs, cex.main = 0.85)
      p <- grDevices::recordPlot()
      p
    })

    subplot_volcano.PLOTLY <- shiny::reactive({
      ngs <- inputData()
      shiny::req(ngs)

      comp <- 1
      gs <- 1
      comp <- gs_contrast()
      ngs <- inputData()
      shiny::req(ngs)

      gxmethods <- selected_gxmethods() ## from module-expression
      shiny::req(gxmethods)

      gx.meta <- ngs$gx.meta$meta[[comp]]
      meta.q <- apply(gx.meta$q[, gxmethods, drop = FALSE], 1, max, na.rm = TRUE)
      limma1 <- data.frame(meta.fx = gx.meta$meta.fx, meta.q = meta.q)
      gx.annot <- ngs$genes[rownames(gx.meta), c("gene_name", "gene_title")]
      limma <- cbind(gx.annot, limma1)

      gs <- gset_selected()
      if (is.null(gs)) {
        return(NULL)
      }
      gs <- gs[1]
      gset <- getGSETS(gs)[[1]]
      jj <- match(toupper(gset), toupper(limma$gene_name))
      sel.genes <- setdiff(limma$gene_name[jj], c(NA, "", " "))

      fdr <- 1
      fdr <- as.numeric(gs_fdr())

      fc.genes <- as.character(limma[, grep("^gene$|gene_name", colnames(limma))])
      fx <- limma[, grep("logFC|meta.fx|fc", colnames(limma))[1]]
      qval <- limma[, grep("^q|adj.P.Val|meta.q|qval|padj", colnames(limma))[1]]

      qval <- pmax(qval, 1e-12) ## prevent q=0
      qval[which(is.na(qval))] <- 1
      xlim <- c(-1, 1) * max(abs(fx), na.rm = TRUE)
      ylim <- c(0, 12)
      ylim <- c(0, max(12, 1.1 * max(-log10(qval), na.rm = TRUE)))
      ylim

      lfc <- 0.20
      lfc <- as.numeric(gs_lfc())
      y <- -log10(qval + 1e-20)
      scaled.fx <- scale(fx, center = FALSE)
      scaled.y <- scale(y, center = FALSE)

      impt <- function(g) {
        j <- match(g, fc.genes)
        x1 <- scaled.fx[j]
        y1 <- scaled.y[j]
        x <- sign(x1) * (x1**2 + 0.25 * y1**2)
        names(x) <- g
        x
      }
      lab.genes <- c(
        head(sel.genes[order(impt(sel.genes))], 10),
        head(sel.genes[order(-impt(sel.genes))], 10)
      )

      plotlyVolcano(
        x = fx, y = y, names = fc.genes,
        source = "plot1", marker.type = "scattergl",
        highlight = sel.genes, label = lab.genes,
        group.names = c("group1", "group0"),
        psig = fdr, lfc = lfc,
        xlab = "effect size (log2FC)",
        ylab = "significance (-log10q)",
        marker.size = 4,
        displayModeBar = FALSE,
        showlegend = FALSE
      ) %>%
        plotly::layout(margin = list(b = 60))
    })

    PlotModuleServer(
      "plot",
      func = subplot_volcano.RENDER,
      func2 = subplot_volcano.PLOTLY,
      plotlib = "base",
      plotlib2 = "plotly",
      pdf.width = 5, pdf.height = 5,
      res = c(72, 100),
      add.watermark = watermark
    )
  })
}
