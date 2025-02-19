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
                                       title,
                                       info.text,
                                       info.methods,
                                       info.references,
                                       info.extra_link,
                                       caption,
                                       height,
                                       width) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("pltmod"),
    label = label,
    plotlib = c("plotly", "ggplot"),
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    title = title,
    caption = caption,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height,
    cards = TRUE,
    card_names = c("dynamic", "static")
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
                                           show_pv,
                                           res,
                                           genes_selected,
                                           labeltype = reactive("symbol"),
                                           watermark = FALSE,
                                           pgx) {
  moduleServer(id, function(input, output, session) {
    # reactive function listening for changes in input
    plot_data <- shiny::reactive({
      shiny::req(res())
      shiny::req(length(comp1()) > 0)
      # calculate required inputs for plotting

      comp1 <- comp1()
      fdr <- as.numeric(fdr())
      lfc <- as.numeric(lfc())
      res <- res()

      symbols <- playbase::probe2symbol(
        probes = rownames(res), res, query = "symbol", fill_na = TRUE
      )

      qval <- pmax(res$meta.q, 1e-20)
      pval <- pmax(res$meta.p, 1e-20)
      x <- res$logFC
      y <- -log10(qval + 1e-12)
      ylab <- "Significance (-log10q)"
      if (show_pv()) {
        y <- -log10(pval + 1e-12)
        ylab <- "Significance (-log10p)"
      }

      names <- ifelse(is.na(res$gene_title), rownames(res), res$gene_title)

      label.names <- playbase::probe2symbol(rownames(res), pgx$genes, labeltype(), fill_na = TRUE)

      return(list(
        x = x,
        y = y,
        ylab = ylab,
        symbols = symbols,
        features = rownames(res),
        names = names,
        sel.genes = genes_selected()$sel.genes,
        lab.genes = genes_selected()$lab.genes,
        fdr = fdr,
        lfc = lfc,
        label.names = label.names
      ))
    })

    plotly.RENDER <- function(marker.size = 4, lab.cex = 1) {
      pd <- plot_data()
      shiny::req(pd)

      plt <- playbase::plotlyVolcano(
        x = pd[["x"]],
        y = pd[["y"]],
        names = pd$features,
        label.names = pd[["label.names"]],
        source = "plot1",
        marker.type = "scattergl",
        highlight = pd[["sel.genes"]],
        label = pd[["lab.genes"]],
        label.cex = lab.cex,
        group.names = c("group1", "group0"),
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        xlab = "Effect size (log2FC)",
        ylab = pd[["ylab"]],
        marker.size = marker.size,
        showlegend = FALSE,
        color_up_down = TRUE
      )
      plt
    }

    modal_plotly.RENDER <- function() {
      fig <- plotly.RENDER(marker.size = 8, lab.cex = 1.5) %>%
        plotly::layout(
          font = list(size = 18),
          legend = list(
            font = list(size = 18)
          )
        )
      fig
    }

    base.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)

      names <- pd$features

      playbase::ggVolcano(
        x = pd[["x"]],
        y = pd[["y"]],
        names = names,
        highlight = pd[["sel.genes"]],
        label = pd[["lab.genes"]],
        label.names = pd[["label.names"]],
        label.cex = 4,
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        xlab = "Effect size (log2FC)",
        ylab = pd[["ylab"]],
        marker.size = 1,
        showlegend = FALSE,
        title = NULL
      )
    }

    base.RENDER.modal <- function() {
      pd <- plot_data()
      shiny::req(pd)

      names <- pd$features

      playbase::ggVolcano(
        x = pd[["x"]],
        y = pd[["y"]],
        names = names,
        highlight = pd[["sel.genes"]],
        label = pd[["lab.genes"]],
        label.names = pd[["label.names"]],
        label.cex = 6,
        axis.text.size = 22,
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        xlab = "Effect size (log2FC)",
        ylab = pd[["ylab"]],
        marker.size = 1.8,
        showlegend = FALSE,
        title = NULL
      )
    }

    plot_data_csv <- function() {
      dt <- plot_data()
      df <- data.frame(dt$x, dt$y)
      colnames(df) <- c("x", "y")
      rownames(df) <- make.unique(dt$symbols)
      return(df)
    }

    plot_grid <- list(
      list(plotlib = "plotly", func = plotly.RENDER, func2 = modal_plotly.RENDER, card = 1),
      list(plotlib = "ggplot", func = base.RENDER, func2 = base.RENDER.modal, card = 2)
    )

    lapply(plot_grid, function(x) {
      PlotModuleServer(
        "pltmod",
        plotlib = x$plotlib,
        func = x$func,
        func2 = x$func2,
        csvFunc = plot_data_csv,
        res = c(80, 95), # resolution of plots
        pdf.width = 10,
        pdf.height = 8,
        add.watermark = watermark,
        card = x$card
      )
    })
  }) ## end of moduleServer
}
