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
    plotlib = "plotly",
    info.text = info.text,
    caption = caption,
    options = plot_options,
    download.fmt = c("png", "pdf", "csv"),
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
expression_plot_volcanoMethods_server <- function(id,
                                                  pgx,
                                                  comp, # input$gx_contrast
                                                  features, # input$gx_features
                                                  fdr, # input$gx_fdr
                                                  lfc, # input$gx_lfc
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## reactive function listening for changes in input
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(comp(), features())
      comp <- comp()
      features <- features()

      fdr <- as.numeric(fdr())
      lfc <- as.numeric(lfc())
      sel.genes <- rownames(pgx$X)
      if (features != "<all>") {
        gset <- playdata::getGSETS(features)
        sel.genes <- unique(unlist(gset))
      }

      pd <- list(
        pgx = pgx,
        fdr = fdr,
        lfc = lfc,
        comp = comp,
        sel.genes = sel.genes
      )

      return(pd)
    })

    plotly_plots <- function(cex = 3, yrange = 0.5, n_rows = 2, margin_l = 50, margin_b = 50) {
      pd <- plot_data()
      shiny::req(pd)

      # Input vars
      sel.genes <- pd[["sel.genes"]]
      fdr <- pd[["fdr"]]
      lfc <- pd[["lfc"]]
      ## meta tables
      comp <- pd[["comp"]]
      mx <- pd[["pgx"]]$gx.meta$meta[[comp]]
      fc <- unclass(mx$fc)
      qv <- unclass(mx$q)
      all_genes <- rownames(mx)
      fc.genes <- pd[["pgx"]]$genes[rownames(mx), "gene_name"]
      nplots <- min(24, ncol(qv))
      methods <- colnames(fc)
      sub_plots <- vector("list", length = length(methods))
      names(sub_plots) <- methods
      if (nplots <= 5) {
        n_rows <- n_rows - 1
        }
      shiny::withProgress(message = "computing volcano plots ...", value = 0, {
        for (i in 1:nplots) {
          # Get plot data
          fx <- fc[, i]
          qval <- qv[, i]
          method_i <- methods[i]
          is.sig <- (qval <= fdr & abs(fx) >= lfc)
          sig.genes <- fc.genes[is.sig]
          qval <- -log(qval + 1e-12)
          title_y <- max(qval) - max(qval) *(yrange/10) 
          # Call volcano plot
          sub_plots[[i]] <- playbase::plotlyVolcano(
            x = fx,
            y = qval,
            names = all_genes,
            source = "plot1",
            marker.type = "scattergl",
            highlight = sig.genes,
            group.names = c("group1", "group0"),
            psig = fdr,
            lfc = lfc,
            marker.size = cex,
            showlegend = FALSE
            # Add plot title
          ) %>% plotly::add_annotations(
              text = paste("<b>", method_i, "</b>"),
              font = list(size = 15),
              showarrow = FALSE,
              xanchor = "left",
              yanchor = "bottom",
              x = 0,
              y = title_y
          )  %>%
            shinyHugePlot::plotly_build_light(.)
        }
      })

      # Pass argument scale_per_plot to subplot
      shareY <- shareX <- ifelse(input$scale_per_plot, TRUE, FALSE)

      # Arrange subplots
      suppressWarnings(
      all_plts <- plotly::subplot(sub_plots, nrows = n_rows , margin = 0.01, 
      titleY = FALSE, titleX = FALSE, shareX = shareX, shareY = shareY) %>%
      
      # Add common axis titles
      plotly::layout(annotations = list(
                list(x = -0.025, y = 0.5, text = "significance (-log10q)",
                     font = list(size = 13),
                     textangle = 270,
                     showarrow = FALSE, xref='paper', yref='paper'),
                list(x = 0.5, y = -0.10, text = "effect size (log2FC)",
                     font = list(size = 13),
                     showarrow = FALSE, xref='paper', yref='paper')
                ),
                  margin = list(l = margin_l, b = margin_b)))

      return(all_plts)
    }

    modal_plotly.RENDER <- function() {
      fig <- plotly_plots(cex = 0.45, yrange = 0.5, n_rows = 2, margin_b = 30)
      return(fig)
    }

    big_plotly.RENDER <- function() {
      fig <- plotly_plots(cex = 0.45, yrange = 0.2, n_rows = 3, margin_b = 85) %>%
        plotly::style(
          marker.size = 6
        )
      return(fig)
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = modal_plotly.RENDER,
      func2 = big_plotly.RENDER,
      res = c(80, 90), ## resolution of plots
      pdf.width = 12, pdf.height = 5,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
