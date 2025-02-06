##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

contrast_correlation_ui <- function(
    id,
    title,
    info.text,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  ctcorrplot.opts <- shiny::tagList(
    ## "Show correlation values in cells."),
    withTooltip(
      shiny::checkboxInput(ns("allfc"), "show all contrasts", FALSE),
      "Show all contrasts or just the selected ones."
    ),
    ##       "Fix heatmap layout when changing number of top genes"),
    withTooltip(
      shiny::radioButtons(ns("ntop"), tspan("number of top genes"),
        c("100", "1000", "all"),
        selected = "1000", inline = TRUE
      ),
      "Number of top genes to compute correlation values."
    )
  )

  PlotModuleUI(
    ns("ctcorrplot"),
    title = title,
    label = "b",
    plotlib = "plotly",
    info.text = info.text,
    caption = caption,
    options = ctcorrplot.opts,
    download.fmt = c("png", "pdf", "csv"),
    height = height,
    width = width
  )
}


contrast_correlation_server <- function(id,
                                        getFoldChangeMatrix,
                                        pgx,
                                        input_comparisons,
                                        watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)

      res <- getFoldChangeMatrix()
      shiny::req(res$fc)
      if (is.null(res)) {
        return(NULL)
      }

      allfc <- input$allfc
      F <- res$fc
      Q <- res$qv

      if (!allfc) {
        F <- F[, input_comparisons(), drop = FALSE]
        Q <- Q[, input_comparisons(), drop = FALSE]
      }

      dbg("[contrast_correlation_server:plot_data] 1: dim(F) = ", dim(F))
      dbg("[contrast_correlation_server:plot_data] 1: dim(Q) = ", dim(Q))

      # module can only run with at least two comparisons
      validate(need(
        dim(F)[2] >= 2,
        "Less than 2 comparisons selected. Please select at least 2 comparison on the settings sidebar."
      ))

      ntop <- 2000
      ntop <- input$ntop
      if (ntop == "all") ntop <- 999999
      ntop <- as.integer(ntop)

      F[is.na(F)] <- 0
      jj <- head(order(-rowMeans(F**2, na.rm = TRUE)), ntop)
      F <- apply(F[jj, , drop = FALSE], 2, rank, na.last = "keep")
      R <- cor(F, use = "pairwise")
      R <- round(R, digits = 2)
      return(R)
    })

    ctcorrplot.PLOTLY <- function() {
      R <- plot_data()
      res <- getFoldChangeMatrix()
      col <- grDevices::colorRampPalette(c(omics_colors("brand_blue"), omics_colors("grey"), omics_colors("red")))(16)
      col <- gplots::colorpanel(64, "royalblue3", "grey90", "indianred3")
      if (min(R, na.rm = TRUE) >= 0) col <- tail(col, 32)
      if (max(R, na.rm = TRUE) <= 0) col <- head(col, 32)

      bluered.pal <- colorRampPalette(colors = c(omics_colors("brand_blue"), omics_colors("grey"), omics_colors("red")))
      cellnote <- NULL

      if (is.null(R) && !is.null(res)) {
        # Deal with ncol == 1
        plt <- plotly::plot_ly()
        plt <- plotly::add_annotations(plt,
          x = 0.5,
          y = 0.5,
          text = "Heatmap plot requires \nmore than two conditions",
          showarrow = FALSE,
          font = list(size = 20)
        )

        plt <- plotly::layout(plt,
          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
        )
      } else {
        plt <- heatmaply::heatmaply(
          R,
          margins = c(250, 200, NA, 0),
          cellnote = cellnote, cellnote_size = 11,
          cellnote_textposition = "middle center",
          colors = bluered.pal,
          limits = c(-1, 1)
        )
      }
      return(plt)
    }

    PlotModuleServer(
      "ctcorrplot",
      func = ctcorrplot.PLOTLY,
      csvFunc = plot_data,
      plotlib = "plotly",
      res = c(80, 85),
      pdf.width = 5, pdf.height = 5,
      add.watermark = watermark
    )
  })
}

# OLD PLOTING FUNCTION
