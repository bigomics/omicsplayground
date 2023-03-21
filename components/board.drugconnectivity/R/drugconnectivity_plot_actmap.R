##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' Activation map plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
drugconnectivity_plot_actmap_ui <- function(id,
                                            label = "",
                                            height = c(750, 1400),
                                            fullH = 750) {
  ns <- shiny::NS(id)
  info_text <- strwrap("The <strong>Activation Matrix</strong> visualizes the
                       activation of drug activation enrichment across the
                       conditions. The size of the circles correspond to their
                       relative activation, and are colored according to their
                       upregulation (red) or downregulation (blue) in the
                       contrast profile.")

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("dsea_normalize"),
        "normalize activation matrix", FALSE
      ),
      "Normalize columns of the activation matrix."
    )
  )
  PlotModuleUI(ns("plot"),
    title = "Activation matrix",
    label = label,
    plotlib = "base",
    info.text = info_text,
    options = plot_opts,
    download.fmt = c("png", "pdf", "csv"),
    height = c(fullH, 750),
    width = c("100%", 1400)
  )
}

#' Activation map plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
drugconnectivity_plot_actmap_server <- function(id,
                                                pgx,
                                                dsea_contrast,
                                                dsea_method,
                                                dsea_table,
                                                getActiveDSEA,
                                                watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      dseaPlotActmap <- function(pgx, dmethod, contr, nterms, nfc) {
        if (is.null(pgx$drugs)) {
          return(NULL)
        }
        nes <- pgx$drugs[[dmethod]]$X
        qv <- pgx$drugs[[dmethod]]$Q

        score <- nes * (1 - qv)**2
        score[is.na(score)] <- 0
        score <- score[order(-score[, contr]**2), , drop = FALSE] ## sort by score

        dsea <- getActiveDSEA()
        res <- dsea$table

        ## filter with table selection/search
        ii <- dsea_table$rows_all()
        shiny::req(ii)
        if (length(ii) > 0) {
          res <- res[ii, , drop = FALSE]
          dd <- intersect(res$drug, rownames(score))
          score <- score[dd, , drop = FALSE]
        }
        if (nrow(score) <= 1) {
          return(NULL)
        }

        score <- head(score, nterms) ## max number of terms
        score <- score[, head(order(-colSums(score**2)), nfc), drop = FALSE] ## max contrs/FC
        score <- score + 1e-3 * matrix(rnorm(length(score)), nrow(score), ncol(score))

        if (input$dsea_normalize) {
          score <- t(t(score) / (1e-8 + sqrt(colMeans(score**2))))
        }
        score <- sign(score) * abs(score)**3 ## fudging
        score <- score / (1e-8 + max(abs(score), na.rm = TRUE))

        if (NCOL(score) > 1) {
          d1 <- as.dist(1 - cor(t(score), use = "pairwise"))
          d2 <- as.dist(1 - cor(score, use = "pairwise"))
          d1[is.na(d1)] <- 1
          d2[is.na(d2)] <- 1
          ii <- hclust(d1)$order
          jj <- hclust(d2)$order
          score <- score[ii, jj, drop = FALSE]
        } else {
          score <- score[order(-score[, 1]), , drop = FALSE]
        }

        colnames(score) <- substring(colnames(score), 1, 30)
        rownames(score) <- substring(rownames(score), 1, 50)
        cex2 <- 0.85
        par(mfrow = c(1, 1), mar = c(1, 1, 1, 1), oma = c(0, 1, 0, 0))

        corrplot::corrplot(score,
          is.corr = FALSE, cl.pos = "n",
          col = BLUERED(100),
          col.lim = c(-1, 1) * max(abs(score), na.rm = TRUE),
          tl.cex = 0.9 * cex2, tl.col = "grey20", tl.srt = 90
        )
      }

      plot_data <- shiny::reactive({
        dsea_contrast <- dsea_contrast()
        dsea_method <- dsea_method()
        shiny::req(pgx, dsea_contrast, dsea_method)
        shiny::validate(shiny::need(
          "drugs" %in% names(pgx),
          "no 'drugs' in object."
        ))
        if (is.null(pgx$drugs)) {
          return(NULL)
        }

        if (is.null(dsea_contrast)) {
          return(NULL)
        }

        res <- list(
          pgx = pgx,
          dsea_contrast = dsea_contrast,
          dsea_method = dsea_method
        )

        return(res)
      })

      plot.RENDER <- shiny::reactive({
        res <- plot_data()
        pgx <- res$pgx
        dsea_contrast <- res$dsea_contrast
        dsea_method <- res$dsea_method

        dseaPlotActmap(pgx, dsea_method, dsea_contrast, nterms = 50, nfc = 20)
      })

      plot.RENDER2 <- shiny::reactive({
        res <- plot_data()
        pgx <- res$pgx
        dsea_contrast <- res$dsea_contrast
        dsea_method <- res$dsea_method

        dseaPlotActmap(pgx, dsea_method, dsea_contrast, nterms = 50, nfc = 100)
      })

      PlotModuleServer(
        "plot",
        plotlib = "base", # does not use plotly
        func = plot.RENDER,
        func2 = plot.RENDER2,
        csvFunc = plot_data,
        res = 72,
        pdf.width = 6, pdf.height = 9,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
