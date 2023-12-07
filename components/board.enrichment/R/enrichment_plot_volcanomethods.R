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
      rm(pd)
      # Call volcano plots
      all_plts <- playbase::plotlyVolcano_multi(FC = fc, 
                                      Q = qv, 
                                      fdr = fdr, 
                                      lfc = lfc,
                                      cex = cex,
                                      title_y =  "significance (-log10q)", 
                                      title_x = "effect size (log2FC)", 
                                      share_axis = !input$scale_per_method,
                                      yrange = yrange,
                                      n_rows = n_rows,
                                      margin_l =  margin_l,
                                      margin_b = margin_b)

      return(all_plts)
    }

    # Render functions
    modal_plotly.RENDER <- function() {
      fig <- plotly_plots(cex = 3, yrange = 0.5, n_rows = 2, margin_b = 20, margin_l = 50) %>%
        playbase::plotly_build_light(.)
      return(fig)
    }

    big_plotly.RENDER <- function() {
      fig <- plotly_plots(yrange = 0.02, n_rows = 3, margin_b = 20, margin_l = 20) %>%
        plotly::style(
          marker.size = 6
        ) %>%
        playbase::plotly_build_light(.)
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
