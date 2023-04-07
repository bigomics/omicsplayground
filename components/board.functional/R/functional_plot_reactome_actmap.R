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
functional_plot_reactome_actmap_ui <- function(
  id,
  title,
  caption,
  info.text,
  label = "",
  height
  ) {
  ns <- shiny::NS(id)
  
  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("reactome_normalize"),
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
    width = c("100%", "100%")
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
functional_plot_reactome_actmap_server <- function(id,
                                               pgx,
                                               getReactomeTable,
                                               watermark = FALSE) {
  moduleServer(
      id, function(input, output, session) {
          
      plotREACTOMEactmap <- function(meta, df,
                                 normalize = 1, nterms = 40, nfc = 10, tl.cex=1) {
        fx <- sapply(meta, function(x) x$meta.fx)
        qv <- sapply(meta, function(x) x$meta.q)
        rownames(fx) <- rownames(qv) <- rownames(meta[[1]])

        kk <- rownames(fx)
        kk <- as.character(df$pathway)

        if (length(kk) < 3) {
          return(NULL)
        }

        if (mean(is.na(qv)) < 0.01) {
          score <- fx[kk, , drop = FALSE] * (1 - qv[kk, , drop = FALSE])**2
        } else {
          score <- fx[kk, , drop = FALSE]
        }
        
        score <- score[head(order(-rowSums(score**2)), nterms), , drop = FALSE] ## nr gene sets
        score <- score[, head(order(-colSums(score**2)), nfc), drop = FALSE] ## max comparisons/FC
        score <- score + 1e-3 * matrix(rnorm(length(score)), nrow(score), ncol(score))
        d1 <- as.dist(1 - cor(t(score), use = "pairwise"))
        d2 <- as.dist(1 - cor(score, use = "pairwise"))
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

        ## fudged score just for visualization
        score2 <- score
        if (normalize) score2 <- t(t(score2) / apply(abs(score2), 2, max))
        score2 <- sign(score2) * abs(score2 / max(abs(score2)))**1 ## fudging
        rownames(score2) <- tolower(gsub(".*:|reactome_|_Homo.*$", "",
          rownames(score2),
          ignore.case = TRUE
        ))
        rownames(score2) <- substring(rownames(score2), 1, 60)
        colnames(score2) <- playbase::shortstring(colnames(score2), 30)
        colnames(score2) <- paste0(colnames(score2), " ")

        bmar <- 0 + pmax(50 - nrow(score2), 0) * 0.3
        par(mfrow = c(1, 1), mar = c(1, 1, 10, 1), oma = c(0, 1.5, 0, 0.5))

        corrplot::corrplot(
          score2,
          is.corr = FALSE,
          cl.pos = "n",
          col = BLUERED(100),
          tl.cex = 1.0*tl.cex,
          tl.col = "grey20",
          tl.srt = 90,
          mar = c(0, 0, 0.5, 0)
        )
      }

      plot_data <- shiny::reactive({
        df <- getReactomeTable()
        res <- list(
          df = df,
          pgx = pgx
        )
        return(res)
      })

      plot_RENDER <- function() {
        res <- plot_data()
        pgx <- res$pgx
        df <- res$df

        if (is.null(df) || nrow(df) == 0) {
          return(NULL)
        }
        meta <- pgx$gset.meta$meta
        plotREACTOMEactmap(
            meta, df,
            normalize = input$reactome_normalize,
            nterms = 25,
            nfc = 25,
            tl.cex = 0.9
        )
      }

      plot_RENDER2 <- function() {
        res <- plot_data()
        pgx <- res$pgx
        df <- res$df

        if (is.null(df) || nrow(df) == 0) {
          return(NULL)
        }
        meta <- pgx$gset.meta$meta
        plotREACTOMEactmap(
            meta, df,
            normalize = input$reactome_normalize,
            nterms = 25,
            nfc = 100,
            tl.cex = 1
        )
      }

      PlotModuleServer(
        "plot",
        plotlib = "base",
        func = plot_RENDER,
        func2 = plot_RENDER2,
        csvFunc = plot_data,
        res = 72,
        pdf.height = 9, pdf.width = 9,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
