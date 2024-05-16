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
    ),
    withTooltip(
      shiny::radioButtons(ns("enrch_volcanoall_subs"), "Subsample plot",
        c("Yes", "No"),
        inline = TRUE, selected = "No"
      ),
      "Number of top genesets to consider for counting the gene frequency."
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

  PlotModuleUI(
    ns("pltmod"),
    plotlib = "plotly",
    title = title,
    info.text = info.text,
    caption = caption,
    options = plot_options,
    height = height,
    width = width,
    download.fmt = c("png", "pdf"),
    subplot = TRUE
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
    shiny::observe({
      shiny::updateSelectInput(session, "pltmod-subplot_selector", choices = subplot_names())
    })

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
      sel.gsets <- rownames(meta[[1]])
      gset_collections <- playbase::pgx.getGeneSetCollections(gsets = rownames(pgx$gsetX))
      sel.gsets <- gset_collections[[gs_features()]]

      # Calc. meta scores and get Q and FC
      FC <- vector("list", length(meta))
      Q <- vector("list", length(meta))
      names(FC) <- names(meta)
      names(Q) <- names(meta)
      for (i in names(meta)) {
        mx <- calcGsetMeta(i, gsmethod, pgx = pgx)
        FC[[i]] <- mx[, "fc", drop = FALSE]
        Q[[i]] <- mx[, "qv", drop = FALSE]
      }

      # Prepare output matrices
      matF <- do.call(cbind, FC)
      matQ <- do.call(cbind, Q)
      colnames(matF) <- names(FC)
      colnames(matQ) <- names(FC)

      pd <- list(
        FC = matF,
        Q = matQ,
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
      fc <- pd$FC
      qv <- pd$Q
      fdr <- pd[["fdr"]]
      lfc <- pd[["lfc"]]
      # Call volcano plots
      all_plts <- playbase::plotlyVolcano_multi(
        FC = fc,
        Q = qv,
        fdr = fdr,
        lfc = lfc,
        cex = cex,
        source = "enrich_volcanoall",
        title_y = "significance (-log10q)",
        title_x = "effect size (log2FC)",
        share_axis = !input$scale_per_method,
        yrange = yrange,
        n_rows = n_rows,
        margin_l = margin_l,
        margin_b = margin_b,
        color_up_down = input$color_up_down
      )


      return(all_plts)
    }

    subplot_names <- shiny::reactive({
      shiny::req(plot_data())
      pd <- plot_data()
      fc <- pd$FC
      names_download <- colnames(fc)
      download_options <- 1:length(colnames(fc))
      names(download_options) <- names_download
      download_options <- c("All", download_options)
      return(download_options)
    })

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
      "pltmod",
      plotlib = "plotly",
      func = modal_plotly.RENDER,
      func2 = big_plotly.RENDER,
      pdf.width = 10,
      pdf.height = 5,
      res = c(72, 85),
      add.watermark = watermark,
      subplot = TRUE
    )
  }) ## end module-server
} ## server
