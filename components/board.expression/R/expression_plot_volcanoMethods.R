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
    title = "Volcano plots for all methods",
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    caption = caption,
    options = plot_options,
    download.fmt = c("png", "pdf", "csv"),
    height = height,
    width = width,
    subplot = TRUE
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
    shiny::observe({
      shiny::updateSelectInput(session, "pltmod-subplot_selector", choices = subplot_names())
    })

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

    plotly_plots <- function(cex = 2, yrange = 0.5, n_rows = 2, margin_l = 50, margin_b = 50) {
      pd <- plot_data()
      shiny::req(pd)

      # Input vars
      sel.genes <- pd[["sel.genes"]]
      fdr <- pd[["fdr"]]
      lfc <- pd[["lfc"]]
      ## meta tables
      comp <- pd[["comp"]]
      mx <- pd[["pgx"]]$gx.meta$meta[[comp]]
      fc <- mx[which(rownames(mx) %in% sel.genes), "fc", drop = FALSE]
      qv <- mx[which(rownames(mx) %in% sel.genes), "q", drop = FALSE]
      rm(mx, pd)
      # Call volcano plots
      all_plts <- playbase::plotlyVolcano_multi(
        FC = fc,
        Q = qv,
        fdr = fdr,
        lfc = lfc,
        cex = cex,
        title_y = "significance (-log10q)",
        title_x = "effect size (log2FC)",
        share_axis = input$scale_per_plot,
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
      sel.genes <- pd[["sel.genes"]]
      comp <- pd[["comp"]]
      mx <- pd[["pgx"]]$gx.meta$meta[[comp]]
      fc <- mx[which(rownames(mx) %in% sel.genes), "fc", drop = FALSE]
      names_download <- colnames(fc[[1]])
      download_options <- 1:length(colnames(fc[[1]]))
      names(download_options) <- names_download
      download_options <- c("All", download_options)
      return(download_options)
    })

    modal_plotly.RENDER <- function() {
      fig <- plotly_plots(cex = 3, yrange = 0.05, n_rows = 2, margin_b = 20, margin_l = 50) %>%
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

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = modal_plotly.RENDER,
      func2 = big_plotly.RENDER,
      csvFunc = plot_data_csv,
      res = c(80, 90), ## resolution of plots
      pdf.width = 12, pdf.height = 5,
      add.watermark = watermark,
      subplot = TRUE
    )
  }) ## end of moduleServer
}
