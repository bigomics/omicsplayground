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
expression_plot_volcano_ui <- function(id,
                                       label = "",
                                       height,
                                       width) {
  ns <- shiny::NS(id)
  options <- tagList(
    actionButton(ns("button1"), "some action")
  )

  info_text <- "<b>Volcano-plot</b> showing fold-change (logFC) versus significance (-log10q) on the x and y axes, respectively. "

  PlotModuleUI(ns("pltmod"),
    title = "Volcano plot",
    label = label,
    plotlib = "plotly",
    ## outputFunc = plotly::plotlyOutput,
    ## outputFunc2 = plotly::plotlyOutput,
    info.text = info_text,
    options = NULL,
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
#' @param comp1
#' @param fdr
#' @param lfc
#' @param features
#' @param res
#' @param sel1
#' @param df1
#' @param watermark
#'
#'
#'
#' @export
expression_plot_volcano_server <- function(id,
                                           comp1,
                                           fdr,
                                           lfc,
                                           features,
                                           res,
                                           sel1,
                                           df1,
                                           sel2,
                                           df2,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    # reactive function listening for changes in input
    plot_data <- shiny::reactive({
      # calculate required inputs for plotting

      comp1 <- comp1()
      fdr <- as.numeric(fdr())
      lfc <- as.numeric(lfc())
      features <- features()
      res <- res()
      sel1 <- sel1()
      df1 <- df1()
      sel2 <- sel2()
      df2 <- df2()

      ## if no gene selected we should show full volcano plot
#      req(sel1())

      fam.genes <- res$gene_name

      if (is.null(res)) {
        return(NULL)
      }
      if (length(comp1) == 0) {
        return(NULL)
      }
      if (is.null(features)) {
        return(NULL)
      }
      if (features != "<all>") {
        gset <- getGSETS(features)
        fam.genes <- unique(unlist(gset))
      }

      jj <- match(toupper(fam.genes), toupper(res$gene_name))
      sel.genes <- res$gene_name[setdiff(jj, NA)]

      fc.genes <- as.character(res[, grep("^gene$|gene_name", colnames(res))])
      qval <- res[, grep("adj.P.Val|meta.q|qval|padj", colnames(res))[1]]
      qval <- pmax(qval, 1e-20)
      x <- res[, grep("logFC|meta.fx|fc", colnames(res))[1]]
      y <- -log10(qval + 1e-12)

      sig.genes <- fc.genes[which(qval <= fdr & abs(x) > lfc)]
      sel.genes <- intersect(sig.genes, sel.genes)
      scaled.x <- scale(x, center = FALSE)
      scaled.y <- scale(y, center = FALSE)
      impt <- function(g) {
        j <- match(g, fc.genes)
        x1 <- scaled.x[j]
        y1 <- scaled.y[j]
        x <- sign(x1) * (x1**2 + 0.25 * y1**2)
        names(x) <- g
        x
      }

      lab.cex <- 1
      gene.selected <- !is.null(sel1) && !is.null(df1)
      gset.selected <- !is.null(sel2) && !is.null(df2)
      if (gene.selected && !gset.selected) {
        lab.genes <- rownames(df1)[sel1]
        sel.genes <- lab.genes
        lab.cex <- 1.3
      } else if (gene.selected && gset.selected) {
        gs <- rownames(df2)[sel2]
        gset <- unlist(getGSETS(gs))
        sel.genes <- intersect(sel.genes, gset)
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
      xlim <- c(-1, 1) * max(abs(x), na.rm = TRUE)
      ylim <- c(0, max(12, 1.1 * max(-log10(qval), na.rm = TRUE)))

      return(list(
        x = x,
        y = y,
        fc.genes = fc.genes,
        sel.genes = sel.genes,
        lab.genes = lab.genes,
        lab.cex = lab.cex,
        fdr = fdr,
        lfc = lfc
      ))
    })


    plotly.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)

      plt <- playbase::plotlyVolcano(
        x = pd[["x"]],
        y = pd[["y"]],
        names = pd[["fc.genes"]],
        source = "plot1",
        marker.type = "scattergl",
        highlight = pd[["sel.genes"]],
        label = pd[["lab.genes"]],
        label.cex = pd[["lab.cex"]],
        group.names = c("group1", "group0"),
        ## xlim=xlim, ylim=ylim, ## hi.col="#222222",
        ## use.fdr=TRUE,
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        xlab = "effect size (log2FC)",
        ylab = "significance (-log10q)",
        marker.size = 3,
        displayModeBar = FALSE,
        showlegend = FALSE
      ) ## %>% plotly::layout(margin = list(b = 35))
      plt
    }

    modal_plotly.RENDER <- function() {    
      fig <- plotly.RENDER() %>%
        plotly::layout(
          font = list(size = 18),
          legend = list(
            font = list(size = 18)
          )
        ) %>% plotly::style(
          marker.size = 6
        )
      fig
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(80, 95), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
