##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Expression plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
expression_plot_maplot_ui <- function(
    id,
    title,
    info.text,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)
  options <- tagList(
    shiny::checkboxInput(
      inputId = ns("color_up_down"),
      label = "Color up/down regulated",
      value = TRUE
    )
  )

  PlotModuleUI(ns("pltmod"),
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    caption = caption,
    options = options,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height
  )
}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#' @param pgx
#' @param gx_fdr
#' @param gx_contrast
#' @param gx_lfc
#' @param gx_features
#' @param res
#' @param sel1
#' @param df1
#' @param watermark
#'
#'
#'
#' @export
expression_plot_maplot_server <- function(id,
                                          pgx,
                                          gx_fdr,
                                          gx_contrast,
                                          gx_lfc,
                                          gx_features,
                                          res,
                                          sel1,
                                          df1,
                                          sel2,
                                          df2,
                                          fam.genes,
                                          genes_in_sel_geneset,
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    # #calculate required inputs for plotting ---------------------------------

    plot_data <- shiny::reactive({
      comp1 <- gx_contrast()
      if (length(comp1) == 0) {
        return(NULL)
      }
      shiny::req(pgx$X)
      sel1 <- sel1()
      df1 <- df1()
      sel2 <- sel2()
      df2 <- df2()
      X <- pgx$X
      res <- res()
      lfc <- as.numeric(gx_lfc())
      fdr <- as.numeric(gx_fdr())

      dbg("[expression_plot_maplot.R] sel1 = ", sel1())

      if (is.null(res)) {
        return(NULL)
      }

      ## filter genes by gene family or gene set
      fam.genes <- res$symbol
      if (gx_features() != "<all>") {
        gset <- playdata::getGSETS(gx_features())
        fam.genes <- unique(unlist(gset))
      }

      jj <- match(toupper(fam.genes), toupper(res$symbol))
      sel.genes <- res$symbol[setdiff(jj, NA)]

      qval <- res[, grep("adj.P.Val|meta.q|qval|padj", colnames(res))[1]]
      y <- res[, grep("logFC|meta.fx|fc", colnames(res))[1]]

      scaled.x <- scale(-log10(qval), center = FALSE)
      scaled.y <- scale(y, center = FALSE)
      fc.genes <- res$symbol

      impt <- function(g) {
        j <- match(g, fc.genes)
        x1 <- scaled.x[j]
        y1 <- scaled.y[j]
        x <- sign(x1) * (0.25 * x1**2 + y1**2)
        names(x) <- g
        x
      }

      sig.genes <- fc.genes[which(qval <= fdr & abs(y) > lfc)]
      sel.genes <- intersect(sig.genes, sel.genes)
      lab.cex <- 1
      gene.selected <- !is.null(sel1) && !is.null(df1)
      gset.selected <- !is.null(sel2) && !is.null(df2)
      if (gene.selected && !gset.selected) {
        lab.genes <- df1$symbol[sel1]
        sel.genes <- lab.genes
        lab.cex <- 1.3
      } else if (gene.selected && gset.selected) {
        sel.genes <- genes_in_sel_geneset()
        lab.genes <- c(
          head(sel.genes[order(impt(sel.genes))], 10),
          head(sel.genes[order(-impt(sel.genes))], 10)
        )
        lab.cex <- 1
      } else {
        lab.genes <- c(
          head(sel.genes[order(impt(sel.genes))], 10),
          head(sel.genes[order(-impt(sel.genes))], 10)
        )
        lab.cex <- 1
      }

      ylim <- c(-1, 1) * max(abs(y), na.rm = TRUE)
      x <- rowMeans(X[rownames(res), ], na.rm = TRUE)

      impt <- function(g) {
        j <- match(g, fc.genes)
        x1 <- scale(x, center = FALSE)[j]
        y1 <- scale(y, center = FALSE)[j]
        x <- sign(y1) * (1.0 * x1**2 + 1.0 * y1**2)
        names(x) <- g
        x
      }
      lab.genes <- c(
        head(sel.genes[order(impt(sel.genes))], 10),
        head(sel.genes[order(-impt(sel.genes))], 10)
      )

      return(list(
        x = x,
        y = y,
        ylim = ylim,
        sel.genes = sel.genes,
        lab.genes = lab.genes,
        lab.cex = lab.cex,
        fc.genes = fc.genes,
        sel1 = sel1,
        df1 = df1,
        sel2 = sel2,
        df2 = df2,
        fdr = fdr,
        lfc = lfc
      ))
    })

    plotly.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)

      par(mfrow = c(1, 1), mar = c(4, 3, 1, 1.5), mgp = c(2, 0.8, 0), oma = c(0, 0, 0, 0))

      plt <- playbase::plotlyMA(
        x = pd[["x"]], y = pd[["y"]], names = pd[["fc.genes"]],
        source = "plot1", marker.type = "scattergl",
        highlight = pd[["sel.genes"]],
        label = pd[["lab.genes"]], label.cex = pd[["lab.cex"]],
        group.names = c("group1", "group0"),
        psig = pd[["fdr"]], lfc = pd[["lfc"]],
        xlab = "average expression (log2.CPM)",
        ylab = "effect size (log2.FC)",
        marker.size = 4,
        displayModeBar = FALSE,
        showlegend = FALSE,
        color_up_down = input$color_up_down
      )
      plt
    }

    modal_plotly.RENDER <- function() {
      fig <- plotly.RENDER() %>%
        plotly::layout(
          font = list(size = 18),
          legend = list(
            font = list(size = 18)
          )
        )
      fig <- plotly::style(fig, marker.size = 10)
      fig
    }

    plot_data_csv <- function() {
      dt <- plot_data()
      df <- data.frame(dt$x, dt$y)
      colnames(df) <- c("x", "y")
      return(df)
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      remove_margins = FALSE,
      csvFunc = plot_data_csv, ##  *** downloadable data as CSV
      res = c(80, 95), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
