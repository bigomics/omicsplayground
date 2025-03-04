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
    plotlib = c("plotly", "ggplot"),
    download.fmt = c("png", "pdf", "svg"),
    cards = TRUE,
    card_names = c("dynamic", "static")
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
    plot_data <- shiny::reactive({
      par(mfrow = c(1, 1), mgp = c(1.2, 0.4, 0), oma = c(0, 0, 0, 0.4))
      par(mar = subplot.MAR)

      shiny::req(pgx$X)
      shiny::validate(shiny::need(!is.null(gset_selected()), tspan("Please select a geneset.", js = FALSE)))

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

      gset <- playbase::probe2symbol(gset, pgx$genes, "gene_name")

      jj <- match(gset, limma$gene_name)
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

      return(list(
        x = fx,
        y = -log10(qval),
        fc.genes = fc.genes,
        sel.genes = sel.genes,
        lab.cex = 1,
        fdr = fdr,
        lfc = lfc
      ))
    })

    plotly.RENDER <- function(marker.size = 3, lab.cex = 1) {
      pd <- plot_data()
      shiny::req(pd)

      playbase::plotlyVolcano(
        x = pd[["x"]],
        y = pd[["y"]],
        names = pd[["fc.genes"]],
        label.names = pd[["fc.genes"]],
        source = "plot1",
        marker.type = "scattergl",
        highlight = pd[["sel.genes"]],
        label = pd[["sel.genes"]],
        label.cex = lab.cex,
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        xlab = "Effect size (log2FC)",
        ylab = "Significance (-log10q)",
        marker.size = marker.size,
        displayModeBar = FALSE,
        showlegend = FALSE,
        color_up_down = TRUE
      ) %>%
        plotly::layout(
          margin = list(l = 0, r = 0, t = 0, b = 0)
        )
    }

    plotly.RENDER2 <- function() {
      plotly.RENDER(marker.size = 8, lab.cex = 1.5) %>%
        plotly::layout(
          font = list(size = 18),
          legend = list(
            font = list(size = 18)
          )
        )
    }

    base.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)

      playbase::ggVolcano(
        x = pd[["x"]],
        y = pd[["y"]],
        title = NULL,
        names = pd[["fc.genes"]],
        label.names = pd[["fc.genes"]],
        highlight = pd[["sel.genes"]],
        label = pd[["sel.genes"]],
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        xlab = "Effect size (log2FC)",
        ylab = "Significance (-log10q)",
        marker.size = 1,
        showlegend = FALSE
      )
    }

    base.RENDER.modal <- function() {
      pd <- plot_data()
      shiny::req(pd)

      playbase::ggVolcano(
        x = pd[["x"]],
        y = pd[["y"]],
        names = pd[["fc.genes"]],
        label.names = pd[["fc.genes"]],
        highlight = pd[["sel.genes"]],
        label = pd[["sel.genes"]],
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        xlab = "Effect size (log2FC)",
        ylab = "Significance (-log10q)",
        marker.size = 2,
        label.cex = 6,
        axis.text.size = 24,
        showlegend = FALSE,
        title = NULL
      )
    }

    plot_grid <- list(
      list(plotlib = "plotly", func = plotly.RENDER, func2 = plotly.RENDER2, card = 1),
      list(plotlib = "ggplot", func = base.RENDER, func2 = base.RENDER.modal, card = 2)
    )

    lapply(plot_grid, function(x) {
      PlotModuleServer(
        "plot",
        plotlib = x$plotlib,
        func = x$func,
        func2 = x$func2,
        res = c(80, 95), # resolution of plots
        pdf.width = 10,
        pdf.height = 8,
        add.watermark = watermark,
        card = x$card
      )
    })
  })
}
