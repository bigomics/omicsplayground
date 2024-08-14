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
#'
#' @export
signature_plot_volcano_ui <- function(
    id,
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
        inputId = ns("share_axis"),
        label = "Share X/Y axis",
        value = TRUE
      ),
      "Share X/Y axis between plots so they show the same range. Easier to compare.",
      placement = "left", options = list(container = "body")
    ),
    withTooltip(
      shiny::checkboxInput(
        inputId = ns("color_up_down"),
        label = "Color up/down regulated",
        value = TRUE
      ),
      "Color up/down regulated features.",
      placement = "left", options = list(container = "body")
    )
  )

  PlotModuleUI(ns("plot"),
    title = title,
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    download.fmt = c("png", "pdf"),
    plotlib = "plotly",
    options = plot_opts,
    height = height,
    width = width
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
signature_plot_volcano_server <- function(id,
                                          pgx,
                                          sigCalculateGSEA,
                                          enrichmentContrastTable,
                                          selected_gxmethods,
                                          enrichmentGeneTable,
                                          getEnrichmentGeneTable,
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ##
    plotly_plots <- function(cex = 3,
                             yrange = 0.5,
                             n_rows = 2,
                             margin_l = 70,
                             margin_b = 50) {
      shiny::req(sigCalculateGSEA())

      ## filter with table selection/search
      ii <- enrichmentContrastTable$rows_selected()
      if (is.null(ii)) {
        ii <- enrichmentContrastTable$rows_all()
      }
      shiny::req(ii)

      # Input vars
      gsea <- sigCalculateGSEA()
      ct <- rownames(gsea$output)[ii]
      mm <- selected_gxmethods()
      meta <- playbase::pgx.getMetaMatrix(pgx, methods = mm)
      fc <- meta$fc[, ct, drop = FALSE]
      qv <- meta$qv[, ct, drop = FALSE]

      features <- rownames(fc)
      symbols <- pgx$genes[rownames(fc), "symbol"]

      # Get gene selected labels
      if (length(ct) == 1) {
        sel.gene <- getEnrichmentGeneTable()$symbol
        row <- enrichmentGeneTable$rows_selected()
        if (length(row)) {
          sel.gene <- sel.gene[row]
        }
      } else {
        sel.gene <- gsea$gset
      }

      share_axis <- input$share_axis

      # Call volcano plots
      all_plts <- playbase::plotlyVolcano_multi(
        FC = fc,
        Q = qv,
        names = features,
        label.names = symbols,
        cex = cex,
        by_sig = FALSE,
        highlight = gsea$gset,
        label = sel.gene,
        share_axis = share_axis,
        yrange = yrange,
        n_rows = n_rows,
        margin_l = margin_l,
        margin_b = margin_b,
        interplot_margin = c(0.02, 0.02, 0.04, 0.04),
        # Remove default titles
        title_y = "",
        title_x = "",
        color_up_down = input$color_up_down
      ) %>%
        plotly::layout(
          annotations = list(
            list(
              x = -0.0,
              xshift = -26,
              xanchor = "right",
              y = 0.5,
              text = "significance (-log10q)",
              font = list(size = 16),
              textangle = 270,
              showarrow = FALSE,
              xref = "paper",
              yref = "paper"
            ),
            list(
              x = 0.5,
              y = 0.0,
              yshift = -20,
              yanchor = "top",
              text = "effect size (log2FC)",
              font = list(size = 16),
              showarrow = FALSE,
              xref = "paper",
              yref = "paper"
            )
          )
        )

      return(all_plts)
    }

    big_plotly.RENDER <- function() {
      nr <- length(enrichmentContrastTable$rows_all())
      n_rows <- floor(sqrt(nr))
      fig <- plotly_plots(
        cex = 5,
        yrange = 0.02,
        n_rows = n_rows,
        margin_b = 45,
        margin_l = 70
      ) %>%
        plotly::style(
          marker.size = 6
        ) %>%
        playbase::plotly_build_light(.)
      return(fig)
    }

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = big_plotly.RENDER,
      res = c(90, 130), ## resolution of plots
      pdf.width = 10,
      pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
