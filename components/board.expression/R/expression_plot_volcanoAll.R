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
expression_plot_volcanoAll_ui <- function(id,
                                          title,
                                          caption,
                                          info.text,
                                          label = "",
                                          height,
                                          width) {
  ns <- shiny::NS(id)

  plot_options <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("scale_per_plot"), "scale per plot", FALSE),
      "Scale each volcano plots individually.",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(
    id = ns("pltmod"),
    title = title,
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
expression_plot_volcanoAll_server <- function(id,
                                              pgx,
                                              getAllContrasts,
                                              features,
                                              fdr,
                                              lfc,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## reactive function listening for changes in input
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(features())

      # Input variables
      features <- features()
      ct <- getAllContrasts()
      FC <- ct$F
      Q <- ct$Q
      fdr <- as.numeric(fdr())
      lfc <- as.numeric(lfc())
      comp <- names(FC)
      shiny::req(length(comp) > 0)

      # Subset genes if requrired
      sel.genes <- rownames(pgx$X)
      if (features != "<all>") {
        gset <- playdata::getGSETS(features)
        sel.genes <- unique(unlist(gset))
      }

      ## combined matrix for output
      matF <- do.call(cbind, FC)
      colnames(matF) <- paste0("fc.", names(FC))
      matQ <- do.call(cbind, Q)
      colnames(matQ) <- paste0("q.", names(Q))
      FQ <- cbind(matF, matQ)

      pd <- list(
        FQ = FQ, ## Remember: the first element is returned as downloadable CSV
        comp = comp,
        fdr = fdr,
        lfc = lfc,
        sel.genes = sel.genes,
        FC = FC,
        Q = Q
      )

      return(pd)
    })

    plotly_plots <- function(cex = 2, yrange = 0.5, n_rows = 2, margin_l = 50, margin_b = 50) {
      pd <- plot_data()
      shiny::req(pd)

      # Input vars
      fdr <- pd[["fdr"]]
      lfc <- pd[["lfc"]]
      ## meta tables
      fc_cols <- grep("fc.*", colnames(pd[["FQ"]]))
      q_cols <- grep("q.*", colnames(pd[["FQ"]]))
      fc <- pd[["FQ"]][, fc_cols, drop = FALSE]
      qv <- pd[["FQ"]][, q_cols, drop = FALSE]
      colnames(fc) <- gsub("fc.", "", colnames(fc))
      colnames(qv) <- gsub("q.", "", colnames(qv))
      rm(pd)
      # Call volcano plots
      all_plts <- playbase::plotlyVolcano_multi(
        FC = fc,
        Q = qv,
        fdr = fdr,
        lfc = lfc,
        cex = cex,
        share_axis = input$scale_per_plot,
        yrange = yrange,
        n_rows = n_rows,
        margin_l = margin_l,
        margin_b = margin_b
      )

      return(all_plts)
    }

    modal_plotly.RENDER <- function() {
      fig <- plotly_plots(cex = 3, yrange = 0.05, n_rows = 2, margin_b = 20, margin_l = 50) %>%
        playbase::plotly_build_light(.)
      return(fig)
    }

    big_plotly.RENDER <- function() {
      fig <- plotly_plots(yrange = 0.2, n_rows = 3, margin_b = 85) %>%
        plotly::style(
          marker.size = 6
        ) %>%
        playbase::plotly_build_light(.)

      return(fig)
    }

    #    shiny::observeEvent(plotly::event_data("plotly_relayout"),{
    #      shiny::req(modal_plotly.RENDER())
    #      ns <- session$ns
    #      print(ns("test"))
    #      shinyHugePlot::updatePlotlyH(session, "expression-volcanoAll-pltmod-renderfigure", plotly::event_data("plotly_relayout"), ds)
    #    })

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = modal_plotly.RENDER,
      func2 = big_plotly.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(70, 90), ## resolution of plots
      pdf.width = 12, pdf.height = 5,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
