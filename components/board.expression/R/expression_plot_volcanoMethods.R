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
expression_plot_volcanoMethods_ui <- function(
    id,
    title,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  plot_options <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("scale_per_plot"), "scale plots", FALSE),
      "Scale each volcano plot individually.",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    title = "Volcano plots for all methods",
    label = label,
    plotlib = c("plotly", "ggplot"),
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    options = plot_options,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    width = width,
    cards = TRUE,
    card_names = c("dynamic", "static")
  )
}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
expression_plot_volcanoMethods_server <- function(id,
                                                  pgx,
                                                  comp, # input$gx_contrast
                                                  fdr, # input$gx_fdr
                                                  lfc, # input$gx_lfc
                                                  show_pv,
                                                  genes_selected,
                                                  labeltype = reactive("symbol"),
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## reactive function listening for changes in input
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(comp())
      comp <- comp()

      fdr <- as.numeric(fdr())
      lfc <- as.numeric(lfc())

      mx <- pgx$gx.meta$meta[[comp]]
      mx.features <- rownames(mx)
      mx.symbols <- pgx$genes[mx.features, "symbol"]

      if (labeltype() == "symbol") {
        names <- mx.features
        label.names <- mx.symbols
      } else {
        names <- mx.symbols
        label.names <- mx.features
      }

      pd <- list(
        mx = mx,
        pgx = pgx,
        fdr = fdr,
        lfc = lfc,
        comp = comp,
        sel.genes = genes_selected()$sel.genes,
        lab.genes = genes_selected()$lab.genes,
        names = names,
        label.names = label.names
      )

      return(pd)
    })

    plotly_plots <- function(cex = 2, yrange = 0.5, n_rows = 2, margin_l = 50, margin_b = 50) {
      pd <- plot_data()
      shiny::req(pd)

      # Input vars
      sel.genes <- pd[["sel.genes"]]
      lab.genes <- pd[["lab.genes"]]
      fdr <- pd[["fdr"]]
      lfc <- pd[["lfc"]]
      names <- pd[["names"]]
      label.names <- pd[["label.names"]]

      ## meta tables
      comp <- pd[["comp"]]
      mx <- pd[["pgx"]]$gx.meta$meta[[comp]]
      x <- mx[, "fc", drop = FALSE]
      y <- mx[, "q", drop = FALSE]
      title_y <- "Significance (-log10q)"
      if (show_pv()) {
        y <- mx[, "p", drop = FALSE]
        title_y <- "Significance (-log10p)"
      }

      mx.features <- rownames(mx)
      mx.symbols <- pgx$genes[mx.features, "symbol"]
      mx.names <- ifelse(is.na(pgx$genes[mx.features, "gene_title"]),
        mx.features,
        pgx$genes[mx.features, "gene_title"]
      )

      label.names <- playbase::probe2symbol(rownames(mx), pgx$genes, labeltype(), fill_na = TRUE)

      # Call volcano plots
      all_plts <- playbase::plotlyVolcano_multi(
        FC = x,
        Q = y,
        fdr = fdr,
        lfc = lfc,
        cex = cex,
        names = mx.features,
        label.names = label.names,
        highlight = sel.genes,
        label = lab.genes,
        title_y = title_y,
        title_x = "Effect size (log2FC)",
        share_axis = !input$scale_per_plot,
        yrange = yrange,
        n_rows = n_rows,
        margin_l = margin_l,
        margin_b = margin_b,
        color_up_down = TRUE,
        by_sig = FALSE
      )
      return(all_plts)
    }


    modal_plotly.RENDER <- function() {
      fig <- plotly_plots(
        cex = 3, yrange = 0.05, n_rows = 2, margin_b = 50, margin_l = 70
      ) %>%
        playbase::plotly_build_light(.)
      return(fig)
    }

    big_plotly.RENDER <- function() {
      fig <- plotly_plots(
        yrange = 0.02, n_rows = 3, margin_b = 70, margin_l = 70
      ) %>%
        plotly::style(
          marker.size = 6
        ) %>%
        playbase::plotly_build_light(.)
      return(fig)
    }

    base.plots <- function(label.cex = 4) {
      pd <- plot_data()
      shiny::req(pd)

      sel.genes <- pd[["sel.genes"]]
      lab.genes <- pd[["lab.genes"]]
      fdr <- pd[["fdr"]]
      lfc <- pd[["lfc"]]
      names <- pd[["names"]]
      label.names <- pd[["label.names"]]

      mx <- pd[["mx"]]
      fc <- mx[, "fc", drop = FALSE][[1]] |> unclass()
      qv <- mx[, "q", drop = FALSE][[1]] |> unclass()
      rownames(fc) <- names
      rownames(qv) <- names

      gene_names <- rep(rownames(fc), each = ncol(fc))
      label.names <- rep(label.names, each = ncol(fc))
      fc <- data.frame(fc) %>%
        tidyr::pivot_longer(
          cols = everything(),
          names_to = "facet",
          values_to = "fc"
        )
      qv <- data.frame(qv) %>%
        tidyr::pivot_longer(
          cols = everything(),
          names_to = "facet",
          values_to = "qv"
        )
      facet <- fc$facet
      x <- fc$fc
      y <- qv$qv
      y <- -log10(y + 1e-12)

      playbase::ggVolcano(
        x = x,
        y = y,
        psig = fdr,
        lfc = lfc,
        facet = facet,
        names = gene_names,
        label.names = label.names,
        highlight = sel.genes,
        label = lab.genes,
        marker.size = 1.2,
        label.cex = label.cex,
        ylab = "Significance (-log10q)",
        xlab = "Effect size (log2FC)",
        showlegend = FALSE,
        title = NULL
      )
    }

    big_base.plots <- function(label.cex = 4) {
      base.plots(label.cex = 5)
    }

    plot_grid <- list(
      list(plotlib = "plotly", func = modal_plotly.RENDER, func2 = modal_plotly.RENDER, card = 1),
      list(plotlib = "ggplot", func = base.plots, func2 = big_base.plots, card = 2)
    )

    plot_data_csv <- function() {
      pd <- plot_data()
      sel.genes <- pd[["sel.genes"]]
      comp <- pd[["comp"]]
      mx <- pd[["pgx"]]$gx.meta$meta[[comp]]
      fc <- mx[which(rownames(mx) %in% sel.genes), "fc", drop = FALSE]
      qv <- mx[which(rownames(mx) %in% sel.genes), "q", drop = FALSE]
      df <- cbind(fc, qv, mx)
      return(df)
    }

    lapply(plot_grid, function(x) {
      PlotModuleServer(
        "pltmod",
        plotlib = x$plotlib,
        func = x$func,
        func2 = x$func2,
        csvFunc = plot_data,
        res = c(80, 90), # resolution of plots
        pdf.width = 12,
        pdf.height = 5,
        add.watermark = watermark,
        card = x$card
      )
    })
  }) ## end of moduleServer
}
