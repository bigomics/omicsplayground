##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

compare_plot_expression_ui <- function(id, label = "", height = c(600, 800)) {
  ns <- shiny::NS(id)

  info_text <- "<b>Multi barplots.</b> Barplots of expression values for multiple comparisons in the two datasets (blue and green). "

  PlotModuleUI(
    ns("plot"),
    title = "Expression",
    label = "a",
    info.text = info_text,
    height = c(440, 700),
    width = c("auto", 1280),
    download.fmt = c("png", "pdf")
  )
}

compare_plot_expression_server <- function(id,
                                           pgx,
                                           dataset2,
                                           input.contrast1,
                                           input.contrast2,
                                           hilightgenes,
                                           getOmicsScoreTable,
                                           score_table,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({

    })

    multibarplot.RENDER <- shiny::reactive({
      pgx1 <- pgx
      pgx2 <- dataset2()

      ct1 <- head(names(pgx1$gx.meta$meta), 3)
      ct2 <- head(names(pgx2$gx.meta$meta), 3)
      ct1 <- input.contrast1()
      ct2 <- input.contrast2()
      shiny::req(ct1)
      shiny::req(ct2)
      if (!all(ct1 %in% names(pgx1$gx.meta$meta))) {
        return(NULL)
      }
      if (!all(ct2 %in% names(pgx2$gx.meta$meta))) {
        return(NULL)
      }

      gene <- "ESR1"
      gene <- "ERBB2"
      genes <- c("ERBB2", "ESR1", "FOXA1", "PGR", "ERBB2", "ESR1", "FOXA1", "PGR")
      genes <- rownames(pgx1$X)
      genes <- hilightgenes()

      df <- getOmicsScoreTable()
      if (is.null(df)) {
        return(NULL)
      }
      sel <- score_table$rows_all() ## from module
      shiny::req(sel)
      genes <- rownames(df)[sel]

      xgenes <- intersect(rownames(pgx1$X), rownames(pgx2$X))
      genes <- head(intersect(genes, xgenes), 8)

      par(mfrow = c(2, 4), mar = c(6, 4, 1, 0), mgp = c(2.0, 0.7, 0), oma = c(0, 0, 0, 0))
      for (gene in genes) {
        x1 <- pgx1$X[gene, ]
        x2 <- pgx2$X[gene, ]
        e1 <- contrastAsLabels(pgx1$model.parameters$exp.matrix[, ct1, drop = FALSE])
        e2 <- contrastAsLabels(pgx2$model.parameters$exp.matrix[, ct2, drop = FALSE])
        m1 <- lapply(e1, function(y) tapply(x1, y, mean))
        m2 <- lapply(e2, function(y) tapply(x2, y, mean))

        grp1 <- paste0("1:", sub(":.*", "", names(m1)))
        grp2 <- paste0("2:", sub(":.*", "", names(m2)))
        grp.names <- c(grp1, grp2)

        b1 <- as.vector(sapply(m1, names))
        b2 <- as.vector(sapply(m2, names))
        bar.names <- c(b1, b2)
        bar.names <- toupper(substring(bar.names, 1, 1))

        srt <- 10
        srt <- 0
        srt <- ifelse(max(nchar(grp.names)) <= 5, 0, 25)

        mm <- cbind(do.call(cbind, m1), do.call(cbind, m2))
        mm.group <- c(rep(1, length(m1)), rep(2, length(m2)))

        gx.barplot(mm,
          srt = srt, main = gene, cex.main = 1.0,
          group = mm.group, cex.names = 0.85,
          group.names = grp.names,
          bar.names = bar.names, voff = 3.5,
          legend = FALSE, cex.legend = 0.9,
          ylab = "expression"
        )
      }
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = multibarplot.RENDER,
      # csvFunc = plot_data,
      pdf.width = 5, pdf.height = 5,
      res = c(95, 130),
      add.watermark = watermark
    )
  })
}
