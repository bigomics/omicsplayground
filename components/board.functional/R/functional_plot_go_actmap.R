##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
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
functional_plot_go_actmap_ui <- function(id,
                                         label = "",
                                         rowH = 660) {
  ns <- shiny::NS(id)
  info_text <- strwrap(
    "The <b>GO activation matrix</b> visualizes the activation of GO terms
    across conditions. From this figure, you can easily detect GO terms that
    are consistently up/down across conditions. The size of the circles
    correspond to their relative activation, and are colored according to their
    upregulation (red) or downregulation (blue) in the contrast profile."
  )

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("go_normalize"),
        "normalize activation matrix",
        FALSE
      ),
      "Click to normalize the columns of the activation matrices."
    )
  )

  PlotModuleUI(ns("plot"),
    title = "Activation matrix",
    label = label,
    plotlib = "base",
    info.text = info_text,
    options = plot_opts,
    height = c(rowH, 750),
    width = c("100%", 1400),
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
      plotGOactmap <- function(score, go, normalize, maxterm, maxfc) {
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

        corrplot::corrplot(score,
          is.corr = FALSE, cl.pos = "n", col = BLUERED(100),
          tl.cex = 0.85, tl.col = "grey20", tl.srt = 90,
          mar = c(bmar, 0, 0, 0)
        )
      }

      plot_data <- shiny::reactive({
        shiny::req(pgx)

        res <- list(
          pgx = pgx
        )
        return(res)
      })

      plot_RENDER <- shiny::reactive({
        res <- plot_data()
        pgx <- res$pgx

        if (is.null(pgx$meta.go)) {
          return(NULL)
        }

        score <- pgx$meta.go$pathscore
        go <- pgx$meta.go$graph

        plotGOactmap(
          score = score, go = go,
          normalize = input$go_normalize,
          maxterm = 50,
          maxfc = 25
        )
      })

      plot_RENDER2 <- shiny::reactive({
        res <- plot_data()
        pgx <- res$pgx

        if (is.null(pgx$meta.go)) {
          return(NULL)
        }

        score <- pgx$meta.go$pathscore
        go <- pgx$meta.go$graph

        plotGOactmap(
          score = score, go = go,
          normalize = input$go_normalize,
          maxterm = 50,
          maxfc = 100
        )
      })

      PlotModuleServer(
        "plot",
        plotlib = "base",
        func = plot_RENDER,
        func2 = plot_RENDER2,
        csvFunc = plot_data,
        res = 72,
        pdf.width = 9, pdf.height = 9,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
