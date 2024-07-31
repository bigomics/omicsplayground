##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

compare_plot_gene_corr_ui <- function(
    id,
    title,
    info.text,
    label = "",
    height = c(600, 800)) {
  ns <- shiny::NS(id)

  genecorr.opts <- shiny::tagList(
    withTooltip(shiny::selectInput(ns("colorby"), "Color by:", choices = NULL, multiple = FALSE),
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
    download.fmt = c("png", "pdf")
  )
}

compare_plot_gene_corr_server <- function(id,
                                          pgx,
                                          dataset2,
                                          input.contrast1,
                                          input.contrast2,
                                          hilightgenes,
                                          getOmicsScoreTable,
                                          score_table,
                                          compute,
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    contrast1 <- shiny::reactiveVal(FALSE)
    contrast2 <- shiny::reactiveVal(FALSE)
    shiny::observeEvent(compute(), {
      contrast1(input.contrast1())
      contrast2(input.contrast2())
    })

    shiny::observeEvent(contrast1(), {
      shiny::req(contrast1())
      ct <- contrast1()
      shiny::updateSelectInput(session, "colorby", choices = ct, selected = ct[1])
    })

    genecorr.RENDER <- shiny::reactive({
      shiny::req(getOmicsScoreTable())
      shiny::req(input$colorby)
      shiny::req(hilightgenes())
      shiny::req(contrast1())
      shiny::req(contrast2())
      shiny::req(score_table())

      pgx1 <- pgx
      pgx2 <- dataset2()
      ct1 <- contrast1()
      ct2 <- contrast2()
      gg <- intersect(rownames(pgx1$X), rownames(pgx2$X))
      kk <- intersect(colnames(pgx1$X), colnames(pgx2$X))

      shiny::validate(shiny::need(
        length(kk) > 0,
        "No common samples between datasets, need at least 10 samples to compute gene correlations."
      ))


      ## conform matrices
      X1 <- pgx1$X[gg, kk]
      X2 <- pgx2$X[gg, kk]
      Y1 <- pgx1$samples[kk, ]
      Y2 <- pgx2$samples[kk, ]
      dset1 <- paste0("1: expression")
      dset2 <- paste0("2: expression")

      df <- getOmicsScoreTable()


      sel <- score_table() ## from module
      shiny::req(sel)

      higenes <- head(rownames(df)[sel], 16)
      shiny::validate(shiny::need(
        length(higenes) > 0,
        "Not valid option."
      ))

      ## Set color for points
      klrpal <- rep(RColorBrewer::brewer.pal(12, "Paired"), 99)
      colorby <- input$colorby
      grp <- playbase::pgx.getContrastGroups(pgx1, colorby, as.factor = TRUE)
      grp <- grp[colnames(X1)]
      klr1 <- klrpal[as.integer(grp)]

      nc <- ceiling(sqrt(length(higenes)))
      nr <- (length(higenes) - 1) %/% nc
      par(mfrow = c(nc, nc), mar = c(2.6, 2.3, 1.5, 0.5) * 1, oma = c(3, 3, 0, 0), mgp = c(2.4, 0.7, 0))
      i <- 1
      for (i in 1:length(higenes)) {
        j <- match(higenes[i], rownames(X1))
        if (is.na(j) | length(j) == 0) {
          X1 <- playbase::rename_by(pgx1$X, pgx1$genes, "symbol")
          X2 <- playbase::rename_by(pgx2$X, pgx2$genes, "symbol")
          j <- match(higenes[i], rownames(X1))
        }
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

    plotly_genecorr.RENDER <- shiny::reactive({
      # Requirements
      shiny::req(getOmicsScoreTable())
      shiny::req(input$colorby)
      shiny::req(hilightgenes())
      shiny::req(contrast1())
      shiny::req(contrast2())
      shiny::req(score_table())

      # Input variables
      pgx1 <- pgx
      pgx2 <- dataset2()
      gg <- intersect(rownames(pgx1$X), rownames(pgx2$X))
      kk <- intersect(colnames(pgx1$X), colnames(pgx2$X))
      dset1 <- paste0("1: expression")
      dset2 <- paste0("2: expression")

      shiny::validate(shiny::need(
        length(kk) > 0,
        "No common samples between datasets, need at least 10 samples to compute gene correlations."
      ))

      ## conform matrices
      X1 <- pgx1$X[gg, kk]
      X2 <- pgx2$X[gg, kk]
      # Get Omics Score Table for high expr. genes
      df <- getOmicsScoreTable()
      sel <- score_table() ## from module
      shiny::req(sel)

      higenes <- head(rownames(df)[sel], 16)
      shiny::validate(shiny::need(
        length(higenes) > 0,
        "Not valid option."
      ))

      ## Set color for points
      klrpal <- rep(1:7, 99)
      klrpal <- rep(RColorBrewer::brewer.pal(12, "Paired"), 99)
      colorby <- input$colorby
      grp <- playbase::pgx.getContrastGroups(pgx1, colorby, as.factor = TRUE)
      grp <- grp[colnames(X1)]
      klr1 <- klrpal[as.integer(grp)]

      # Aseemble subplots
      sub_plots <- vector("list", length(higenes))
      names(sub_plots) <- higenes
      for (gene_i in higenes) {
        # Get genes and titles
        j <- match(gene_i, rownames(X1))
        title_y <- max(X2[j, ]) + max(X2[j, ]) * .1
        show_legend <- gene_i == higenes[1]
        xtitle <- ifelse(gene_i %in% higenes[15], dset1, "")
        ytitle <- ifelse(gene_i %in% higenes[9], dset2, "")
        plt <- plotly::plot_ly() %>%
          # Axis
          plotly::layout(
            xaxis = list(title = xtitle, font = list(size = 5)),
            yaxis = list(title = ytitle, font = list(size = 5)),
            legend = list(x = 0.5, y = -0.1, xanchor = "center", orientation = "h", bgcolor = "transparent")
          ) %>%
          # Add the points
          plotly::add_trace(
            x = X1[j, ], y = X2[j, ], name = grp, color = klr1, type = "scatter", mode = "markers",
            marker = list(size = 10),
            showlegend = show_legend
          ) %>%
          plotly::add_annotations(
            text = paste("<b>", gene_i, "</b>"),
            font = list(size = 10),
            showarrow = FALSE,
            xanchor = "left",
            yanchor = "bottom",
            x = 1.5,
            y = title_y
          )
        sub_plots[[gene_i]] <- plt
      }

      # Assemble all subplot in to grid
      suppressWarnings(
        plt <- plotly::subplot(sub_plots,
          nrows = 4, margin = 0.03,
          titleX = TRUE, titleY = TRUE
        )
      )
      return(plt)
    })

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = plotly_genecorr.RENDER, # genecorr.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(80, 90),
      add.watermark = watermark
    )
  })
}
