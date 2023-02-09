##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
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
expression_plot_topfoldchange_ui <- function(id,
                                             label = "",
                                             height,
                                             width) {
  ns <- shiny::NS(id)
  info_text <- "The fold change summary barplot across all contrasts for a gene that is selected from the differential expression analysis table under the <code>Table</code> section."

  PlotModuleUI(ns("pltmod"),
    title = "Gene in contrasts",
    label = label,
    plotlib = "plotly",
    info.text = info_text,
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
#' @param comp
#' @param ngs
#' @param sel
#' @param res
#' @param watermark
#'
#'
#'
#' @export
expression_plot_topfoldchange_server <- function(id,
                                                 comp,
                                                 ngs,
                                                 sel,
                                                 res,
                                                 watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    # #calculate required inputs for plotting ---------------------------------

    plot_data <- shiny::reactive({
      comp <- comp() # input$gx_contrast
      ngs <- ngs()
      sel <- sel()
      res <- res()

      psel <- rownames(res)[sel]
      gene <- ngs$genes[psel, "gene_name"]

      if (is.null(sel) || length(sel) == 0) { # Ugly
        return(list(sel = sel))
      }

      if (is.null(comp) || length(comp) == 0) {
        return(NULL)
      }
      fc <- sapply(ngs$gx.meta$meta, function(x) x[psel, "meta.fx"])
      top.up <- head(names(sort(fc[which(fc > 0)], decreasing = TRUE)), 10)
      top.dn <- head(names(sort(fc[which(fc < 0)], decreasing = FALSE)), 10)
      fc.top <- c(fc[top.up], fc[top.dn])
      fc.top <- fc.top[head(order(-abs(fc.top)), 15)]
      fc.top <- sort(fc.top)
      fc.top <- head(c(fc.top, rep(NA, 99)), 15)

      klr.pal <- RColorBrewer::brewer.pal(4, "Paired")[2:1]
      klr <- klr.pal[1 + 1 * (sign(fc.top) < 0)]

      return(list(
        sel = sel,
        fc.top = fc.top,
        klr = klr,
        gene = gene
      ))
    })

    plotly.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)

      if (is.null(pd[["sel"]]) || length(pd[["sel"]]) == 0) {
        frame()
        text(0.5, 0.5, "No gene selected", col = "black")
        return(NULL)
      }

      if (is.null(res) || is.null(sel)) {
        return(NULL)
      }

      par(mfrow = c(1, 1), mar = c(4, 4, 2, 2) * 1, mgp = c(2, 0.8, 0), oma = c(1, 1, 1, 0.5) * 0.2)
      par(mfrow = c(1, 1), mar = c(6, 3, 0, 1), mgp = c(2, 0.8, 0), oma = c(1, 0, 0, 0))
      nch <- max(nchar(names(pd[["fc.top"]])))
      m1 <- ifelse(nch > 12, 12, 8)
      m1 <- ifelse(nch > 30, 16, m1)

      par(mar = c(3.2, m1 - 0.5, 1, 1))
      cex1 <- 0.9
      nn <- sum(!is.na(pd[["fc.top"]]))
      if (nn > 15) cex1 <- 0.8



      # barplot(pd[["fc.top"]],
      #   col = pd[["klr"]], horiz = TRUE, las = 1,
      #   xlim = c(-1, 1) * max(abs(pd[["fc.top"]]), na.rm = TRUE),
      #   cex.names = cex1, xlab = "fold change (log2)"
      # )
      # title(pd[["gene"]], cex.main = 1, line = -0.15)

      pgx.barplot.PLOTLY(
        data = data.frame(
          x = names(pd[["fc.top"]]),
          y = as.numeric(pd[["fc.top"]])
        ),
        x = "x",
        y = "y",
        title = pd[["gene"]],
        yaxistitle = "Fold change (log2)",
        xaxistitle = "Groups",
        yrange = c(-1.1, 1.1) * max(abs(pd[["fc.top"]])),
        fillcolor = pd[["klr"]],
        plotRawValues = TRUE
      )
    }



    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      # func2 = modal_plotly.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(80, 95), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
