##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' Single cell plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
singlecell_plot_crosstabPlot_ui <- function(id,
                                            label = "",
                                            height,
                                            width,
                                            parent) {
  ns <- shiny::NS(id)

  crosstab.opts <- shiny::tagList(
    withTooltip(shiny::selectInput(parent("crosstabvar"), label = "x-axis:", choices = NULL, multiple = FALSE),
      "Choose a predefined phenotype group on the x-axis.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(shiny::selectInput(parent("crosstabpheno"), label = "y-axis:", choices = NULL, multiple = FALSE),
      "Choose a predefined phenotype group on the y-axis.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(shiny::selectInput(parent("crosstabgene"), label = "gene:", choices = NULL, multiple = FALSE),
      "Visualize the expression barplot of a gene by specifying the gene name.",
      placement = "top", options = list(container = "body")
    )
    ## checkboxGroupInput('crosstaboptions','',c("gene"), inline=TRUE, width='50px')
    ## selectInput("crosstabgene",label=NULL, choices=NULL, multiple=FALSE),
    ## br(), cellArgs=list(width='80px')
  )

  crosstabModule_info <- "The <strong>Proportions tab</strong> visualizes the interrelationships between two categorical variables (so-called cross tabulation). Although this feature is very suitable for a single-cell sequencing data, it provides useful information about the proportion of different cell types in samples obtained by the bulk sequencing method."


  PlotModuleUI(ns("plot"),
    label = label,
    info.text = crosstabModule_info,
    title = "Proportions",
    options = crosstab.opts,
    download.fmt = c("png", "pdf", "csv"),
    height = height,
    width = width
  )
}

#' Single cell plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @export
singlecell_plot_crosstabPlot_server <- function(id,
                                                inputData,
                                                samplefilter,
                                                getDeconvResults2,
                                                crosstabvar, # input$crosstabvar
                                                pheno, # input$crosstabpheno
                                                gene, # input$crosstabgene
                                                watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- shiny::reactive({
      ## if(!input$tsne.all) return(NULL)

      ngs <- inputData()
      crosstabvar <- crosstabvar()
      gene <- gene()
      pheno <- pheno()

      scores <- ngs$deconv[[1]][[1]] ## just an example...
      if (crosstabvar == "<cell type>") {
        scores <- getDeconvResults2()
        if (is.null(scores)) {
          return(NULL)
        }
        scores <- pmax(scores, 0) ## ??
      } else {
        x <- as.character(ngs$Y[, 1])
        x <- as.character(ngs$Y[, crosstabvar])
        x[is.na(x)] <- "_"
        scores <- model.matrix(~ 0 + x)
        rownames(scores) <- rownames(ngs$Y)
        colnames(scores) <- sub("^x", "", colnames(scores))
      }

      ## restrict to selected sample set
      kk <- head(1:nrow(scores), 1000)
      kk <- 1:nrow(scores)
      kk <- selectSamplesFromSelectedLevels(ngs$Y, samplefilter())
      scores <- scores[kk, , drop = FALSE]
      scores <- scores[, which(colSums(scores) > 0), drop = FALSE]
      scores[which(is.na(scores))] <- 0
      dim(scores)

      ## limit to top25??
      topsel <- head(order(-colSums(scores)), 25)
      scores <- scores[, topsel]

      message(
        "[SingleCellBoard::crosstab.plotFUNC] 2 : dim(scores) = ",
        paste(dim(scores), collapse = "x"), "\n"
      )

      message(
        "[SingleCellBoard::crosstab.plotFUNC] 2 : length(kk) = ",
        length(kk), "\n"
      )

      ## expected counts per stat level
      ## kk.counts <- colSums(ngs$counts[,kk,drop=FALSE])  ## total count of selected samples
      kk.counts <- colSums(2**ngs$X[, kk, drop = FALSE]) ## approximate counts from log2X
      grp.counts <- (t(scores / rowSums(scores)) %*% matrix(kk.counts, ncol = 1))[, 1]

      getProportionsTable <- function(pheno, is.gene = FALSE) {
        y <- NULL
        ## if("gene" %in% input$crosstaboptions) {
        if (is.gene) {
          xgene <- ngs$genes[rownames(ngs$X), ]$gene_name
          gx <- ngs$X[which(xgene == pheno), kk, drop = FALSE]
          gx.highTH <- mean(gx, na.rm = TRUE)
          y <- paste(pheno, c("low", "high"))[1 + 1 * (gx >= gx.highTH)]
          table(y)
        } else if (pheno %in% colnames(ngs$samples)) {
          y <- ngs$samples[kk, 1]
          y <- ngs$samples[kk, pheno]
          pheno <- tolower(pheno)
        } else if (pheno == "<cell type>") {
          res1 <- getDeconvResults2()
          res1 <- pmax(res1, 0) ## ??
          res1 <- res1[kk, , drop = FALSE]
          ## res1 <- res1[,which(colSums(res1)>0),drop=FALSE]
          y <- colnames(res1)[max.col(res1)] ## take maximum col??
          remove(res1)
          pheno <- "<cell type>"
        } else {
          return(NULL)
        }

        ## calculate proportions by group
        grp <- factor(as.character(y))
        ngrp <- length(levels(grp))
        grp.score <- apply(scores, 2, function(x) tapply(x, grp, mean, na.rm = TRUE))
        ngrp
        if (ngrp == 1) {
          grp.score <- matrix(grp.score, nrow = 1)
          rownames(grp.score) <- y[1]
          colnames(grp.score) <- colnames(scores)
        }

        ## weighted counts
        grp.score[is.na(grp.score)] <- 0
        grps <- levels(grp)
        grp.score
        fy <- (table(y) / sum(!is.na(y)))
        jj <- match(rownames(grp.score), names(fy))
        grp.score <- grp.score * as.vector(fy[jj])
        ## normalize to total 100%
        grp.score <- grp.score / (1e-20 + sum(grp.score))
        dim(grp.score)

        ## reduce to maximum number of items (x-axis)
        if (0 && ncol(grp.score) > 25) {
          ## jj <- which(colSums(grp.score) > 0.001)
          jj <- order(-colSums(grp.score))
          j1 <- head(jj, 25) ## define maximum number of items
          j0 <- setdiff(jj, j1)
          grp.score0 <- grp.score[, j1, drop = FALSE]
          grp.counts0 <- grp.counts[j1]
          if (length(j0) > 0) {
            grp.score0 <- cbind(grp.score0, "other" = rowSums(grp.score[, j0, drop = FALSE]))
            grp.counts0 <- c(grp.counts0, "other" = sum(grp.counts[j0]))
          }
          grp.score <- grp.score0
          grp.counts <- grp.counts0
          grp.score <- t(t(grp.score) / (1e-20 + colSums(grp.score))) ##
          dim(grp.score)
        }

        ## normalize to total 100% and reduce to maximum number of items (y-axis)
        if (nrow(grp.score) > 10) {
          jj <- order(-rowSums(grp.score))
          j1 <- head(jj, 10) ## define maximum number of items
          j0 <- setdiff(jj, j1)
          grp.score0 <- grp.score[j1, , drop = FALSE]
          if (length(j0) > 0) {
            grp.score0 <- rbind(grp.score0, "other" = colSums(grp.score[j0, , drop = FALSE]))
          }
          grp.score <- grp.score0
        }
        grp.score <- t(t(grp.score) / (1e-20 + colSums(grp.score))) ##


        ## cluster columns??
        ## dist1 <- dist(t(scale(grp.score)))
        dist1 <- dist(t(grp.score))
        dist1[is.na(dist1)] <- mean(dist1, na.rm = TRUE)
        jj <- hclust(dist1)$order
        grp.score <- grp.score[, jj, drop = FALSE]
        return(grp.score)
      }

      ## select phenotype variable
      head(ngs$samples)
      pheno <- 1
      pheno <- "cluster"
      pheno <- "activated"
      pheno <- "cell.type"
      pheno <- "<cell type>"
      pheno <- pheno
      if (is.null(pheno)) {
        return(NULL)
      }

      ## pheno="cluster"
      grp.score1 <- getProportionsTable(pheno, is.gene = FALSE)
      grp.score2 <- NULL
      gene <- ngs$genes$gene_name[1]
      gene <- gene
      if (gene != "<none>") {
        grp.score2 <- getProportionsTable(pheno = gene, is.gene = TRUE)
        kk <- colnames(grp.score2)[order(grp.score2[1, ])]
        grp.score2 <- grp.score2[, match(kk, colnames(grp.score2))]
        grp.score1 <- grp.score1[, match(kk, colnames(grp.score1))]
      }

      jj <- match(colnames(grp.score1), names(grp.counts))
      grp.counts <- grp.counts[jj] / 1e6 ## per million
      names(grp.counts) <- colnames(grp.score1)

      return(list(
        grp.score2 = grp.score2,
        grp.counts = grp.counts,
        grp.score1 = grp.score1
      ))
    })

    plot.render <- function() {
      pd <- plot_data()

      grp.score2 <- pd[["grp.score2"]]
      grp.counts <- pd[["grp.counts"]]
      grp.score1 <- pd[["grp.score1"]]
      ## -------------- plot by estimated cell.type ----------------------

      ## par(mar = c(4,6,2,3))
      plotly::layout(matrix(c(1, 2, 3), 3, 1), heights = c(2, 4, 3))
      if (!is.null(grp.score2)) {
        plotly::layout(matrix(c(1, 2, 3, 4), 4, 1), heights = c(2.2, 1, 4, 2))
      }

      ## top bar with counts
      par(mar = c(0, 5, 5.3, 3), mgp = c(2.0, 0.8, 0))
      xlim <- c(0, 1.2 * length(grp.counts)) ## reserves space for legend
      barplot(grp.counts,
        col = "grey50", ylab = "counts (M)", cex.axis = 0.8,
        cex.lab = 0.8, names.arg = NA, xpd = NA, xlim = 1.3 * xlim, ## log="y",
        ylim = c(0.01, max(grp.counts)), yaxs = "i"
      )
      ## title(pheno, cex.main=1.2, line=2, col="grey40")

      ## middle plot (gene)
      if (!is.null(grp.score2)) {
        klrpal2 <- COLORS[1:nrow(grp.score2)]
        klrpal2 <- rev(grey.colors(nrow(grp.score2), start = 0.45))
        par(mar = c(0, 5, 0.3, 3), mgp = c(2.4, 0.9, 0))
        barplot(100 * grp.score2,
          col = klrpal2, las = 3, srt = 45,
          xlim = 1.3 * xlim, ylim = c(0, 99.99),
          names.arg = rep(NA, ncol(grp.score2)),
          ylab = "(%)", cex.axis = 0.90
        )
        legend(1.02 * xlim[2], 100,
          legend = rownames(grp.score2),
          fill = klrpal2, xpd = TRUE, cex = 0.8, y.intersp = 0.8,
          bg = "white", bty = "n"
        )
      }

      if (1) {
        ## main proportion graph
        klrpal1 <- COLORS[1:nrow(grp.score1)]
        par(mar = c(4, 5, 0.3, 3), mgp = c(2.4, 0.9, 0))
        barplot(100 * grp.score1,
          col = klrpal1, las = 3, srt = 45, xlim = 1.3 * xlim,
          ylim = c(0, 99.99), ylab = "proportion (%)", cex.axis = 0.90
        )
        legend(1.02 * xlim[2], 100,
          legend = rev(rownames(grp.score1)),
          fill = rev(klrpal1), xpd = TRUE, cex = 0.8, y.intersp = 0.8,
          bg = "white", bty = "n"
        )
      }
    }

    PlotModuleServer(
      "plot",
      func = plot.render,
      res = c(110, 110),
      pdf.width = 12, pdf.height = 8,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
