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
correlation_plot_scattercorr_ui <- function(
  id,
  title,
  info.text,
  caption,
  height,
  width) {
  ns <- shiny::NS(id)

  cor_scatter.opts <- shiny::tagList(
    withTooltip(shiny::selectInput(ns("cor_group"), "Color by:", choices = NULL, multiple = FALSE),
      "Variable to split and color by groups.",
      placement = "top"
    ),
    shiny::checkboxInput(ns("corscatter.swapaxis"), "Swap axes")
  )

  PlotModuleUI(ns("plot"),
    title = title,
    plotlib = "base",
    label = "c",
    info.text = info.text,
    caption = caption,
    options = cor_scatter.opts,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height
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
correlation_plot_scattercorr_server <- function(id,
                                                getFilteredExpression,
                                                pgx,
                                                getPartialCorrelationMatrix,
                                                getGeneCorr,
                                                cor_gene,
                                                COL,
                                                watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    shiny::observe({
      px <- colnames(pgx$Y)
      s1 <- grep("^[.]", px, value = TRUE, invert = TRUE)[1]
      shiny::updateSelectInput(session, "cor_group", choices = px, selected = s1)
    })

    cor_scatter.DATA <- shiny::reactive({
      cor_gene <- cor_gene()
      shiny::req(cor_gene)

      X <- getFilteredExpression()
      this.gene <- rownames(pgx$X)[1]
      this.gene <- cor_gene

      R <- getGeneCorr()
      dbg("[cor_scatter.DATA] 1: dim.R = ", dim(R))
      rho <- R[,"cor"]

      return(list(rho, this.gene))
    })

    cor_scatter.PLOTFUN <- function() {
      dt <- cor_scatter.DATA()
      rho <- dt[[1]]
      this.gene <- dt[[2]]
      if (length(rho) == 0) {
        return(NULL)
      }

      colorby <- input$cor_group
      shiny::req(colorby)

      ph <- factor(pgx$samples[, colorby])
      klrpal <- rep(COL, 99)
      klr <- klrpal[as.integer(ph)]

      ndim <- ncol(pgx$X)

      cex_levels <- c(1.2, 0.8, 0.5, 0.2)
      dim_cuts <- c(0, 40, 100, 200, Inf)
      cex <- cex_levels[findInterval(ndim, dim_cuts)]

      par(
        mfrow = c(5, 5), mar = c(4, 3.5, 0.3, 0),
        mgp = c(1.9, 0.7, 0), oma = c(0, 0, 0.5, 0.5)
      )
      par(
        mfrow = c(5, 5), mar = c(3, 3.5, 0.5, 0),
        mgp = c(1.7, 0.6, 0), oma = c(0, 0, 0.5, 0.5)
      )

      swapaxis <- FALSE
      swapaxis <- TRUE
      swapaxis <- input$corscatter.swapaxis
      xlab <- this.gene

      i <- 1
      for (i in 1:min(25, length(rho))) {
        gene2 <- names(rho)[i]
        if (swapaxis) {
          x <- pgx$X[gene2, ]
          y <- pgx$X[this.gene, ]
          xlab <- gene2
          ylab <- this.gene
        } else {
          y <- pgx$X[gene2, ]
          x <- pgx$X[this.gene, ]
          ylab <- gene2
          xlab <- this.gene
        }
        base::plot(x, y,
          pch = 19, cex = cex, col = klr,
          ylab = ylab, xlab = xlab
        )

        y <- y + 1e-3 * rnorm(length(y))
        x <- x + 1e-3 * rnorm(length(x))
        abline(lm(y ~ x), col = "black", lty = 2)

        if (i %% 5 == 1) {
          tt <- c("   ", levels(ph))
          legend("topleft",
            legend = tt,
            fill = c(NA, klrpal), inset = c(0.02, 0.02),
            border = c(NA, rep("black", 99)),
            cex = 0.9, box.lwd = 0, pt.lwd = 0,
            x.intersp = 0.5, y.intersp = 0.8
          )
          legend("topleft", colorby,
            x.intersp = -0.2,
            cex = 0.9, y.intersp = 0.45, bty = "n"
          )
        }
      }
      plot <- recordPlot()
      return(plot)
    }

    PlotModuleServer(
      "plot",
      plotlib = "base",
      func = cor_scatter.PLOTFUN,
      csvFunc = cor_scatter.DATA, ##  NOTE: Not sure what should be the plot data!
      res = c(80, 95), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
