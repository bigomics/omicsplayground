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
    height = height,
    width = c("auto", "100%"),
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
    multibarplot.RENDER <- shiny::reactive({
      shiny::req(input.contrast1())
      shiny::req(input.contrast2())
      shiny::req(getOmicsScoreTable())
      shiny::req(hilightgenes())
      shiny::req(score_table())
      pgx1 <- pgx
      pgx2 <- dataset2()

      org1 <- playbase::pgx.getOrganism(pgx1)
      org2 <- playbase::pgx.getOrganism(pgx2)
      if (org1 == "human") {
        org1 <- "Human"
      } else if(org1 == "mouse") {
        org1 <- "Mouse"
      } 
      if (org2 == "human") {
        org2 <- "Human"
      } else if(org2 == "mouse") {
        org2 <- "Mouse"
      } 

      ct1 <- head(names(pgx1$gx.meta$meta), 3)
      ct2 <- head(names(pgx2$gx.meta$meta), 3)
      ct1 <- input.contrast1()
      ct2 <- input.contrast2()
      if (is.null(pgx1$version) && is.null(pgx2$version)) {
        target_col1 <- target_col2<- "gene_name"
        } else if (org1 == org2) {
        target_col1 <- target_col2<- "symbol"
        } else if (org1 != org2) {
        target_col1 <- target_col2<- "human_ortholog"
        if(!target_col1 %in% colnames(pgx1$genes)) target_col1 <- "gene_name"
        if(!target_col2 %in% colnames(pgx2$genes)) target_col2 <- "gene_name"
      }

      X1 <- playbase::rename_by(pgx1$X, pgx1$genes, target_col1)
      X2 <- playbase::rename_by(pgx2$X, pgx2$genes, target_col2)

      genes <- rownames(X1)
      genes <- hilightgenes()

      df <- getOmicsScoreTable()

      sel <- score_table() ## from module
      genes <- rownames(df)[sel]

      xgenes <- intersect(rownames(X1), rownames(X2))
      genes <- head(intersect(genes, xgenes), 8)

      par(mfrow = c(2, 4), mar = c(6, 4, 1, 0), mgp = c(2.0, 0.7, 0), oma = c(0, 0, 0, 0))
      for (gene in genes) {
        x1 <- X1[gene, ]
        x2 <- X2[gene, ]
        e1 <- playbase::contrastAsLabels(pgx1$model.parameters$exp.matrix[, ct1, drop = FALSE])
        e2 <- playbase::contrastAsLabels(pgx2$model.parameters$exp.matrix[, ct2, drop = FALSE])
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
        playbase::gx.barplot(mm,
          srt = srt, main = gene, cex.main = 1.0,
          group = mm.group, cex.names = 0.85,
          group.names = grp.names,
          bar.names = bar.names, voff = 3.5,
          legend = FALSE, cex.legend = 0.9,
          ylab = "expression"
        )
      }
      p <- grDevices::recordPlot()
      return(p)
    })

    PlotModuleServer(
      "plot",
      func = multibarplot.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(95, 130),
      add.watermark = watermark
    )
  })
}
