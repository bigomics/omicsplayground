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
  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        inputId = ns("color_up_down"),
        label = "Color up/down regulated",
        value = TRUE
      ),
      "Color up/down regulated features.",
      placement = "left", options = list(container = "body")
    ),
    withTooltip(
      shiny::checkboxInput(
        inputId = ns("show_p_values"),
        label = "Plot nominal p-values on the y-axis",
        value = FALSE
      ),
      "Plot nominal p-values on the y-axis.",
      placement = "left", options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    options = plot_opts,
    title = title,
    caption = caption,
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
                                           res,
                                           genes_selected,
                                           watermark = FALSE) {
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

      fc.genes <- playbase::probe2symbol(probes = rownames(res), res, query = "symbol", fill_na = TRUE)

      qval <- res[, grep("adj.P.Val|meta.q|qval|padj", colnames(res))[1]]
      qval <- pmax(qval, 1e-20)
      pval <- res[, grep("pvalue|meta.p|pval|p|p_value", colnames(res))[1]]
      pval <- pmax(pval, 1e-20)
      x <- res[, grep("logFC|meta.fx|fc", colnames(res))[1]]
      y <- -log10(qval + 1e-12)
      y.lab <- "Significance (-log10q)"
      if (input$show_p_values) {
        y <- -log10(pval + 1e-12)
        y.lab <- "Significance (-log10p)"
      }

      return(list(
        x = x,
        y = y,
        y.lab = y.lab,
        fc.genes = fc.genes,
        sel.genes = genes_selected()$sel.genes,
        lab.genes = genes_selected()$lab.genes,
        lab.cex = 1,
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
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        xlab = "Effect size (log2FC)",
        ylab = pd[["y.lab"]],
        marker.size = 3,
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
        ) %>%
        plotly::style(
          marker.size = 6
        )
      fig
    }

    plot_data_csv <- function() {
      dt <- plot_data()
      df <- data.frame(dt$x, dt$y)
      colnames(df) <- c("x", "y")
      rownames(df) <- dt$fc.genes
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
