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
  info.methods,
  info.extra_link,
  caption,
  label = "",
  height,
  width
) {
  ns <- shiny::NS(id)

  options <- tagList()

  PlotModuleUI(ns("pltmod"),
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    caption = caption,
    options = NULL,
    outputFunc = plotly::plotlyOutput,
    outputFunc2 = plotly::plotlyOutput,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "scatter_updown"
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
                                          genes_selected,
                                          labeltype = reactive("symbol"),
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      comp1 <- gx_contrast()
      if (length(comp1) == 0) {
        return(NULL)
      }
      shiny::req(pgx$X)

      X <- pgx$X
      res <- res()
      lfc <- as.numeric(gx_lfc())
      fdr <- as.numeric(gx_fdr())

      if (is.null(res)) {
        return(NULL)
      }

      y <- res[, grep("logFC|meta.fx|fc", colnames(res))[1]]
      ylim <- c(-1, 1) * max(abs(y), na.rm = TRUE)
      x <- rowMeans(X[rownames(res), ], na.rm = TRUE)
      symbols <- rownames(res)

      names <- ifelse(is.na(res$gene_title), rownames(res), res$gene_title)
      label.names <- pgx$genes[rownames(res), ]$gene_name

      shape <- rep("circle", nrow(res))
      names(shape) <- rownames(res)
      if (any(is.na(pgx$counts)) && !any(is.na(pgx$X))) {
        jj <- which(!is.na(pgx$contrasts[, comp1]))
        if (any(jj)) {
          counts <- pgx$counts[rownames(res), rownames(pgx$contrasts)[jj], drop = FALSE]
          nas <- apply(counts, 1, function(x) sum(is.na(x)))
          na.features <- names(nas)[which(nas > 0)]
          shape[match(na.features, names(shape))] <- "cross"
        }
      }

      jj <- which(rownames(res) == genes_selected()$sel.genes)
      if (any(jj)) {
        fc <- y[jj]
        shiny::validate(shiny::need(!is.na(fc), "Fold change for this feature is NA"))
      }

      plot_data <- list(
        x = x,
        y = y,
        ylim = ylim,
        sel.genes = genes_selected()$sel.genes,
        lab.genes = genes_selected()$lab.genes,
        symbols = symbols,
        features = rownames(res),
        names = names,
        shape = shape,
        fdr = fdr,
        lfc = lfc,
        label.names = label.names
      )

      return(plot_data)
    })


    plotly.RENDER <- function(marker.size = 4, lab.cex = 1) {
      pd <- plot_data()
      shiny::req(pd)

      names <- pd[["features"]]
      label.names <- pd[["label.names"]]

      col_up   <- if (!is.null(input$color_up))   input$color_up   else get_color_theme()$primary
      col_down <- if (!is.null(input$color_down)) input$color_down else get_color_theme()$secondary

      lab.genes <- pd[["lab.genes"]]
      if (isTRUE(input$custom_labels) && !is.null(input$label_features) && nchar(trimws(input$label_features)) > 0) {
        lab.genes <- trimws(strsplit(input$label_features, "[,\n]+")[[1]])
        lab.genes <- lab.genes[lab.genes != ""]
      }

      highlight <- pd[["sel.genes"]]
      if (isTRUE(input$color_selection) && length(lab.genes) > 0) {
        highlight <- lab.genes
      }

      plt <- playbase::plotlyMA(
        x = pd[["x"]],
        y = pd[["y"]],
        names = names,
        label.names = label.names,
        highlight = highlight,
        label = lab.genes,
        label.cex = lab.cex,
        shape = pd[["shape"]],
        group.names = c("group1", "group0"),
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        xlab = "Average expression (log2)",
        ylab = "Effect size (log2FC)",
        marker.size = marker.size,
        displayModeBar = FALSE,
        showlegend = FALSE,
        source = "plot1",
        marker.type = "scattergl",
        color_up_down = TRUE,
        colors = c(up = col_up, notsig = "#8F8F8F", down = col_down)
      )
      plt
    }

    modal_plotly.RENDER <- function() {
      fig <- plotly.RENDER(marker.size = 8, lab.cex = 1.5) %>%
        plotly::layout(font = list(size = 18), legend = list(font = list(size = 18)))
      fig
    }

    plot_data_csv <- function() {
      dt <- plot_data()
      df <- data.frame(dt$features, dt$symbols, dt$x, dt$y)
      colnames(df) <- c("feature", "symbol", "x", "y")
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
      add.watermark = watermark,
      download.contrast.name = gx_contrast,
      parent_session = session
    )
  }) ## end of moduleServer
}
