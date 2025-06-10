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
    card_names = c("dynamic", "static"),
    editor = TRUE,
    ns_parent = ns,
    plot_type = "volcano"
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
                                           pval_cap,
                                           pgx) {

  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(res())
      shiny::req(length(comp1()) > 0)

      comp1 <- comp1()
      fdr <- as.numeric(fdr())
      lfc <- as.numeric(lfc())
      res <- res()

      symbols <- playbase::probe2symbol(
        probes = rownames(res), res, query = "symbol", fill_na = TRUE
      )

      pval_cap <- pval_cap()
      qval <- pmax(res$meta.q, pval_cap)
      pval <- pmax(res$meta.p, pval_cap)
      x <- res$logFC
      y <- -log10(qval + pval_cap)
      ylab <- "Significance (-log10q)"
      if (show_pv()) {
        y <- -log10(pval + pval_cap)
        ylab <- "Significance (-log10p)"
      }

      names <- ifelse(is.na(res$gene_title), rownames(res), res$gene_title)
      label.names <- playbase::probe2symbol(rownames(res), pgx$genes, labeltype(), fill_na = TRUE)

      shape <- rep("circle", nrow(res))
      names(shape) <- rownames(res)
      if (any(is.na(pgx$counts)) && !any(is.na(pgx$X))) {
        jj <- which(!is.na(pgx$contrasts[, comp1]))
        if (any(jj)) {
          counts <- pgx$counts[rownames(res), rownames(pgx$contrasts)[jj], drop = FALSE]
          nas <- apply(counts, 1, function(x) sum(is.na(x)))
          na.features <- names(nas)[which(nas>0)]
          shape[match(na.features, names(shape))] <- "cross"
        }
      }
      
      return(list(
        x = x,
        y = y,
        ylab = ylab,
        symbols = symbols,
        features = rownames(res),
        names = names,
        shape = shape,
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
        shape = pd[["shape"]],
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
      
      if (input$custom_labels) {
        label_features <- if (input$label_features == "") {
          NULL
        } else {
          strsplit(input$label_features, "\\s+")[[1]]
        }
        label_features <- label_features
      } else {
        label_features <- pd[["lab.genes"]]
      }

      highlight <- if (input$color_selection) {
        label_features
      } else {
        pd[["sel.genes"]]
      }

      p <- playbase::ggVolcano(
        x = pd[["x"]],
        y = pd[["y"]],
        names = names,
        highlight = highlight,
        label = label_features,
        label.names = pd[["label.names"]],
        label.cex = input$label_size,
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        xlab = "Effect size (log2FC)",
        ylab = pd[["ylab"]],
        marker.size = input$marker_size,
        showlegend = FALSE,
        title = NULL,
        axis.text.size = input$axis_text_size,
        colors = c(
          up = input$color_up,
          notsig = "#707070AA",
          notsel = "#cccccc88",
          down = input$color_down
        )
      )

      if (input$margin_checkbox) {
        margin_top <- ifelse(is.na(input$margin_top), 10, input$margin_top)
        margin_right <- ifelse(is.na(input$margin_right), 10, input$margin_right)
        margin_bottom <- ifelse(is.na(input$margin_bottom), 10, input$margin_bottom)
        margin_left <- ifelse(is.na(input$margin_left), 10, input$margin_left)
        p <- p + ggplot2::theme(
          plot.margin = ggplot2::margin(
            t = margin_top,
            r = margin_right,
            b = margin_bottom,
            l = margin_left,
            unit = "pt"
          )
        )
      }
      
      if (input$aspect_ratio_checkbox) {
        if (is.na(input$aspect_ratio)) {
          p <- p + ggplot2::theme(
            aspect.ratio = 0.5
          )
        } else {
          p <- p + ggplot2::theme(
            aspect.ratio = input$aspect_ratio
          )
        }
      }

      p
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
        card = x$card,
        parent_session = session,
        download.contrast.name = comp1
      )
    })
  }) ## end of moduleServer
}
