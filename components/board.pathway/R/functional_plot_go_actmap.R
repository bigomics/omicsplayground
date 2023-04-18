##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Importance plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
functional_plot_go_actmap_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height,
  width
  ) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("normalize"),
        "normalize activation matrix",
        FALSE
      ),
      "Click to normalize the columns of the activation matrices."
    )
  )

  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    caption = caption,
    plotlib = "base",
    info.text = info.text,
    options = plot_opts,
    height = height,
    width = width,
  )
}

#' Importance plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
functional_plot_go_actmap_server <- function(id,
                                             pgx,
                                             watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
          
      plotGOactmap <- function(score, go, normalize, maxterm, maxfc,
                                 tl.cex=0.85)
        {
        rownames(score) <- igraph::V(go)[rownames(score)]$Term

        ## avoid errors!!!
        score[is.na(score) | is.infinite(score)] <- 0
        score[is.na(score)] <- 0

        ## reduce score matrix
        score <- score[head(order(-rowSums(score**2, na.rm = TRUE)), maxterm), , drop = FALSE] ## max number terms
        score <- score[, head(order(-colSums(score**2, na.rm = TRUE)), maxfc), drop = FALSE] ## max comparisons/FC
        score <- score + 1e-3 * matrix(rnorm(length(score)), nrow(score), ncol(score))

        ## normalize colums
        if (normalize) {
          ## column scale???
          score <- t(t(score) / (1e-8 + sqrt(colMeans(score**2, na.rm = TRUE))))
        }
        score <- score / max(abs(score), na.rm = TRUE) ## global normalize
        score <- sign(score) * abs(score)**0.5 ## fudging for better colors

        d1 <- as.dist(1 - cor(t(score), use = "pairwise"))
        d2 <- as.dist(1 - cor(score, use = "pairwise"))
        d1 <- dist(score)
        d2 <- dist(t(score))
        d1[is.na(d1)] <- 1
        d2[is.na(d2)] <- 1
        ii <- 1:nrow(score)
        jj <- 1:ncol(score)
        if (NCOL(score) == 1) {
          score <- score[order(-score[, 1]), 1, drop = FALSE]
        } else {
          ii <- hclust(d1)$order
          jj <- hclust(d2)$order
          score <- score[ii, jj, drop = FALSE]
        }

        colnames(score) <- substring(colnames(score), 1, 30)
        rownames(score) <- substring(rownames(score), 1, 50)
        colnames(score) <- paste0(colnames(score), " ")

        bmar <- 0 + pmax((50 - nrow(score)) * 0.25, 0)
        par(mfrow = c(1, 1), mar = c(1, 1, 1, 1), oma = c(0, 1.5, 0, 0.5))

        corrplot::corrplot(
          score,
          is.corr = FALSE,
          cl.pos = "n",
          col = BLUERED(100),
          tl.cex = tl.cex,
          tl.col = "grey20",
          tl.srt = 90,
          mar = c(bmar, 0, 0, 0)
        )
      }

      plot_data <- shiny::reactive({
        shiny::req(pgx$meta.go)
        pathscore = pgx$meta.go$pathscore
        graph = pgx$meta.go$graph
        res <- list(
            pathscore = pathscore,
            graph = graph
        )
      })

      plot_RENDER <- function() {
        res <- plot_data()
        shiny::req(res)
        pathscore <- res$pathscore
        graph <- res$graph

        plotGOactmap(
          score = pathscore,
          go = graph,
          normalize = input$normalize,
          maxterm = 50,
          maxfc = 25,
          tl.cex = 1.1
        )
      }

      plot_RENDER2 <- function() {
        res <- plot_data()
        shiny::req(res)

        pathscore <- res$pathscore
        graph <- res$graph

        plotGOactmap(
          score = pathscore,
          go = graph,
          normalize = input$normalize,
          maxterm = 50,
          maxfc = 100,
          tl.cex = 0.85          
        )
      }

      PlotModuleServer(
        "plot",
        plotlib = "base",
        func = plot_RENDER,
        func2 = plot_RENDER2,
        csvFunc = plot_data,
        res = 72,
        pdf.width = 9,
        pdf.height = 9,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
