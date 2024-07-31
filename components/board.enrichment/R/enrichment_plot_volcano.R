##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_volcano_ui <- function(
    id,
    title,
    caption,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    height,
    width) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        inputId = ns("color_up_down"),
        label = "Color up/down regulated",
        value = TRUE
      ),
      "Color up/down regulated features.",
      placement = "left", options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = "a",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    height = height,
    width = width,
    caption = caption,
    plotlib = "plotly",
    plotlib2 = "plotly",
    options = plot_opts,
    download.fmt = c("png", "pdf")
  )
}

enrichment_plot_volcano_server <- function(id,
                                           pgx,
                                           gs_contrast,
                                           selected_gxmethods,
                                           gset_selected,
                                           gs_fdr,
                                           gs_lfc,
                                           subplot.MAR,
                                           geneDetails,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    volcano.RENDER <- shiny::reactive({
      par(mfrow = c(1, 1), mgp = c(1.2, 0.4, 0), oma = c(0, 0, 0, 0.4))
      par(mar = subplot.MAR)

      shiny::req(pgx$X)
      shiny::validate(shiny::need(!is.null(gset_selected()), "Please select a geneset."))

      comp <- 1
      gs <- 1
      comp <- gs_contrast()
      shiny::req(pgx$X)

      gxmethods <- selected_gxmethods() ## from module-expression
      shiny::req(gxmethods)

      gx.meta <- pgx$gx.meta$meta[[comp]]
      meta.q <- apply(gx.meta$q[, gxmethods, drop = FALSE], 1, max) ## max q-value
      limma1 <- data.frame(meta.fx = gx.meta$meta.fx, meta.q = meta.q)
      gx.annot <- pgx$genes[rownames(gx.meta), c("gene_name", "gene_title")]
      limma <- cbind(gx.annot, limma1)

      gset <- geneDetails()$feature

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

      lfc <- 0.20
      lfc <- as.numeric(gs_lfc())

      playbase::plotlyVolcano(
        x = fx,
        y = -log10(qval),
        names = fc.genes,
        source = "plot1",
        marker.type = "scattergl",
        highlight = sel.genes,
        label = sel.genes,
        psig = fdr,
        lfc = lfc,
        xlab = "Effect size (log2FC)",
        ylab = "Significance (-log10q)",
        marker.size = 3,
        displayModeBar = FALSE,
        showlegend = FALSE,
        color_up_down = input$color_up_down
      ) %>%
        plotly::layout(margin = list(b = 60))
    })

    PlotModuleServer(
      "plot",
      func = volcano.RENDER,
      func2 = volcano.RENDER,
      plotlib = "plotly",
      plotlib2 = "plotly",
      pdf.width = 5, pdf.height = 5,
      res = c(72, 100),
      add.watermark = watermark
    )
  })
}
