##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_volcanoall_ui <- function(
    id,
    title,
    info.text,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  plot_options <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("scale_per_method"), "scale per method", FALSE),
      "Scale the volcano plots individually per method..",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("plot"),
    plotlib = "plotly",
    title = title,
    info.text = info.text,
    caption = caption,
    options = plot_options,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

enrichment_plot_volcanoall_server <- function(id,
                                              pgx,
                                              gs_features,
                                              gs_statmethod,
                                              gs_fdr,
                                              gs_lfc,
                                              calcGsetMeta,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(gs_features())

      meta <- pgx$gset.meta$meta
      gsmethod <- colnames(meta[[1]]$fc)
      gsmethod <- gs_statmethod()
      if (is.null(gsmethod) || length(gsmethod) == 0) {
        return(NULL)
      }

      fdr <- as.numeric(gs_fdr())
      lfc <- as.numeric(gs_lfc())
      sel.gsets <- NULL
      sel.gsets <- rownames(meta[[1]])
      gset_collections <- playbase::pgx.getGeneSetCollections(gsets = rownames(pgx$gsetX))
      sel.gsets <- gset_collections[[gs_features()]]

      i <- 1
      mx.list <- list()
      for (i in 1:length(meta)) {
        mx.list[[i]] <- calcGsetMeta(i, gsmethod, pgx = pgx)
      }
      names(mx.list) <- names(meta)

      Q <- lapply(mx.list, function(mx) mx[, "qv"])
      FC <- lapply(mx.list, function(mx) mx[, "fc"])
      names(FC) <- names(Q) <- names(mx.list)

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
      FC <- pd$FC
      Q <- pd$Q
      fdr <- pd[["fdr"]]
      lfc <- pd[["lfc"]]
      ## meta tables
      nplots <- min(24, length(pd$Q))
      sub_plots <- vector("list", length = length(nplots))
      if (nplots <= 5) {
        n_rows <- n_rows - 1
        }
      shiny::withProgress(message = "computing volcano plots ...", value = 0, {
        for (i in 1:nplots) {
          # Get plot data
          fx <- FC[[i]]
          qval <- Q[[i]]
          fc.genes <- names(qval)
          cond_i <- names(pd[["Q"]])[i]
          is.sig1 <- (qval <= fdr & abs(fx) >= lfc)
          sig.genes <- names(fx)[which(is.sig1)]
          qval <- -log(qval + 1e-12)
          title_y <- max(qval) - max(qval) *(yrange/10) 
          # Call volcano plot
          sub_plots[[i]] <- playbase::plotlyVolcano(
            x = fx,
            y = qval,
            names = fc.genes,
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
              text = paste("<b>", cond_i, "</b>"),
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
      "plot",
      plotlib = "plotly",
      func = modal_plotly.RENDER,
      func2 = big_plotly.RENDER,
      pdf.width = 10,
      pdf.height = 5,
      res = c(72, 85),
      add.watermark = watermark
    )
  }) ## end module-server
} ## server
