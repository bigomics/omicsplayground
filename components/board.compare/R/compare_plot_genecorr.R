##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

compare_plot_genecorr_ui <- function(
  id,
  title,
  info.text,
  label = "",
  height = c(600, 800)
) {
  ns <- shiny::NS(id)

  genecorr.opts <- shiny::tagList(
    withTooltip(
      shiny::selectInput(
        ns("colorby"), "Color by:",
        choices = NULL, multiple = FALSE
      ),
      "Color samples by phenotype.",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("plot"),
    plotlib = "plotly",
    title = title,
    label = "c",
    info.text = info.text,
    options = genecorr.opts,
    height = height,
    width = c("auto", "100%"),
    download.fmt = c("png", "pdf", "svg")
  )
}

compare_plot_genecorr_server <- function(id,
                                         pgx,
                                         dataset2,
                                         contrast1,
                                         contrast2,
                                         getScoreTable,
                                         getMatrices,
                                         selected,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    shiny::observeEvent(contrast1(), {
      shiny::req(contrast1())
      ct <- contrast1()
      shiny::updateSelectInput(
        session, "colorby",
        choices = ct, selected = ct[1]
      )
    })

    plot.RENDER <- shiny::reactive({
      shiny::req(getScoreTable())
      shiny::req(input$colorby)
      shiny::req(contrast1())
      shiny::req(contrast2())
      shiny::req(selected())

      pgx1 <- pgx
      pgx2 <- dataset2()
      ct1 <- contrast1()
      ct2 <- contrast2()

      mat <- getMatrices()
      X1 <- mat$X1
      X2 <- mat$X2
      kk <- intersect(colnames(X1), colnames(X2))

      shiny::validate(shiny::need(
        length(kk) >= 10,
        tspan("No common samples between datasets, need at least 10 samples to compute gene correlations.", js = FALSE)
      ))

      ## conform matrices
      gg <- intersect(rownames(X1), rownames(X2))
      X1 <- X1[gg, kk]
      X2 <- X2[gg, kk]
      Y1 <- pgx1$samples[kk, ]
      Y2 <- pgx2$samples[kk, ]
      dset1 <- paste0("1: expression")
      dset2 <- paste0("2: expression")

      df <- getScoreTable()
      shiny::req(selected())
      higenes <- head(selected(), 16)

      ## Set color for points
      klrpal <- rep(RColorBrewer::brewer.pal(12, "Paired"), 99)
      colorby <- input$colorby
      grp <- playbase::pgx.getContrastGroups(pgx1, colorby, as.factor = FALSE)
      grp <- grp[colnames(X1)]
      klr1 <- klrpal[as.integer(grp)]

      nc <- ceiling(sqrt(length(higenes)))
      nr <- (length(higenes) - 1) %/% nc
      par(mfrow = c(nc, nc), mar = c(2.6, 2.3, 1.5, 0.5) * 1, oma = c(3, 3, 0, 0), mgp = c(2.4, 0.7, 0))
      i <- 1
      for (i in 1:length(higenes)) {
        j <- match(higenes[i], rownames(X1))

        base::plot(X1[j, ], X2[j, ],
          xlab = "", ylab = "",
          pch = 20, col = klr1, cex = 1.2
        )
        title(higenes[i], line = 0.4, cex.main = 1.1)
        if (i %% nc == 1) {
          mtext(dset2, 2, line = 2, cex = 0.8)
        }
        if ((i - 1) %/% nc == nr) {
          mtext(dset1, 1, line = 2, cex = 0.8)
        }

        if (i %% nc == 1) {
          tt <- c("   ", levels(grp))
          legend("topleft",
            legend = tt,
            fill = c(NA, klrpal),
            border = c(NA, "black", "black"), bty = "n",
            cex = 0.92, box.lwd = 0, pt.lwd = 0,
            x.intersp = 0.5, y.intersp = 0.8
          )
          legend("topleft", colorby,
            x.intersp = -0.2,
            cex = 0.92, y.intersp = 0.45, bty = "n"
          )
        }
      }
      p <- grDevices::recordPlot()
      p
    })

    plotly.RENDER <- shiny::reactive({
      # Requirements
      shiny::req(getScoreTable())
      shiny::req(input$colorby)
      shiny::req(contrast1())
      shiny::req(contrast2())
      shiny::req(selected())

      # Input variables
      pgx1 <- pgx
      pgx2 <- dataset2()
      mat <- getMatrices()
      X1 <- mat$X1
      X2 <- mat$X2
      kk <- intersect(colnames(X1), colnames(X2))

      shiny::validate(shiny::need(
        length(kk) >= 8,
        tspan("No common samples between datasets, need at least 8 samples to compute gene correlation.", js = FALSE)
      ))

      dset1 <- tspan(paste("expression in", pgx1$name), js = FALSE)
      dset2 <- tspan(paste("expression in", pgx2$name), js = FALSE)

      ## conform matrices
      X1 <- X1[, kk]
      X2 <- X2[, kk]

      # Get Omics Score Table for high expr. genes
      df <- getScoreTable()

      shiny::req(selected())
      higenes <- head(selected(), 16)

      ## Set color for points
      klrpal <- rep(1:7, 99)
      klrpal <- rep(RColorBrewer::brewer.pal(12, "Paired"), 99)
      colorby <- input$colorby
      grp <- playbase::pgx.getContrastGroups(pgx1, colorby, as.factor = FALSE)
      grp <- grp[colnames(X1)]
      grp[is.na(grp)] <- "_"
      grp <- factor(grp)
      klr1 <- klrpal[as.integer(grp)]

      # Assemble subplots
      sub_plots <- vector("list", length(higenes))
      names(sub_plots) <- higenes
      for (gene_i in higenes) {
        # Get genes and titles
        j <- match(gene_i, rownames(X1))

        title_x <- min(X1[j, ], na.rm = TRUE)
        title_y <- max(X2[j, ], na.rm = TRUE)
        show_legend <- (gene_i == higenes[1])
        xtitle <- ytitle <- ""
        ##        xtitle <- ifelse(gene_i %in% higenes[15], dset1, "")
        ##        ytitle <- ifelse(gene_i %in% higenes[9], dset2, "")

        plt <- plotly::plot_ly() %>%
          # Add the points
          plotly::add_trace(
            x = X1[j, ],
            y = X2[j, ],
            name = grp,
            text = colnames(X1),
            color = klr1,
            type = "scatter",
            mode = "markers",
            marker = list(size = 6),
            showlegend = show_legend
          ) %>%
          plotly::add_annotations(
            text = paste("<b>", gene_i, "</b>"),
            font = list(size = 12),
            showarrow = FALSE,
            xanchor = "left",
            yanchor = "bottom",
            x = title_x,
            y = title_y,
            yshift = -0
          ) %>%
          # Axis
          plotly::layout(
            xaxis = list(title = xtitle, font = list(size = 5)),
            yaxis = list(title = ytitle, font = list(size = 5)),
            legend = list(
              x = 0.5, y = -0.1, xanchor = "center",
              orientation = "h", bgcolor = "transparent"
            ),
            plot_bgcolor = "#fcfcff"
          )

        sub_plots[[gene_i]] <- plt
      }

      # Assemble all subplot in to grid
      plt <- plotly::subplot(
        sub_plots,
        nrows = min(4, length(sub_plots)),
        margin = 0.03,
        titleX = TRUE,
        titleY = TRUE
      )

      ## add axis title
      plt <- plt %>%
        plotly::layout(
          annotations = list(
            list(
              x = -0.0,
              xshift = -26,
              xanchor = "right",
              y = 0.5,
              text = dset2,
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
              text = dset1,
              font = list(size = 16),
              showarrow = FALSE,
              xref = "paper",
              yref = "paper"
            )
          )
        )
      return(plt)
    })

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = plotly.RENDER, # genecorr.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(80, 90),
      add.watermark = watermark
    )
  })
}
