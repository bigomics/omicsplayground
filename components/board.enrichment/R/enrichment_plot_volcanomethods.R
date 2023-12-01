##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_volcanomethods_ui <- function(
    id,
    title,
    info.text,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)


  plot_options <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("scale_per_method"), "scale per method", TRUE),
      "Scale the volcano plots individually per method..",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("plot"),
    plotlib = "plotly",
    title = title,
    caption = caption,
    options = plot_options,
    info.text = info.text,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

enrichment_plot_volcanomethods_server <- function(id,
                                                  pgx,
                                                  gs_features,
                                                  gs_contrast,
                                                  gs_fdr,
                                                  gs_lfc,
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx$X, gs_features(), gs_contrast())

      cmp <- gs_contrast()
      mx <- pgx$gset.meta$meta[[cmp]]
      FC <- unclass(mx$fc)
      Q <- unclass(mx$q)
      rownames(Q) <- rownames(FC) <- rownames(mx)
      FC[which(is.infinite(FC))] <- NA
      Q[which(is.na(Q))] <- 1

      fdr <- as.numeric(gs_fdr())
      lfc <- as.numeric(gs_lfc())
      gset_collections <- playbase::pgx.getGeneSetCollections(gsets = rownames(pgx$gsetX))
      sel.gsets <- gset_collections[[gs_features()]]
      nlq <- -log10(1e-99 + unlist(Q))

      pd <- list(
        FC = FC,
        Q = Q,
        sel.gsets = sel.gsets,
        fdr = fdr,
        lfc = lfc
      )
      pd
    })

    plotly_plots <- function(cex = 3, yrange = 0.5, n_rows = 2, margin_l = 50, margin_b = 50) {
      pd <- plot_data()
      shiny::req(pd)

      # Input vars
      fdr <- pd[["fdr"]]
      lfc <- pd[["lfc"]]
      ## meta tables
      fc <- pd$FC
      qv <- pd$Q
      sel.gsets <- pd$sel.gsets

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
          sig.genes <- names(fx)[which(is.sig)]
          if (!is.null(sel.gsets)) sig.genes <- intersect(sel.gsets, sig.genes)
          qval <- -log(qval + 1e-12)
          title_y <- max(qval) - max(qval) *(yrange/10) 
          # Call volcano plot
          sub_plots[[i]] <- playbase::plotlyVolcano(
            x = fx,
            y = qval,
            names = names(fx),
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
      shareY <- shareX <- ifelse(input$scale_per_method, TRUE, FALSE)

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
                list(x = 0.5, y = -0.10, text = "difference  (log2FC)",
                     font = list(size = 13),
                     showarrow = FALSE, xref='paper', yref='paper')
                ),
                  margin = list(l = margin_l, b = margin_b)))

      return(all_plts)
    }

    # Render functions
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
      "plot",
      plotlib = "plotly",
      func = modal_plotly.RENDER,
      func2 = big_plotly.RENDER,
      pdf.width = 10, pdf.height = 5,
      res = c(75, 90),
      add.watermark = watermark
    )
  })
}
