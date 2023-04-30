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
    withTooltip(
      shiny::selectInput(
        ns("colorby"), "Color by:", choices = NULL, multiple = FALSE),
      "Variable to split and color by groups.",
      placement = "top"
    ),
    withTooltip(
      shiny::radioButtons(ns("layout"), "Layout:", c("3x3", "4x4", "5x5"),
        selected="4x4", inline = TRUE),
      "Choose the layout for correlation plots.",
    ),
    withTooltip(
      shiny::checkboxInput(ns("swapaxis"), "Swap XY-axes"),
      "Transpose plot, i.e. swap X-axis and Y-axis.",
    )
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
                                                cor_table,
                                                getPartialCorrelationMatrix,
                                                getGeneCorr,
                                                cor_gene,
                                                COL,
                                                watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    shiny::observe({
      px <- colnames(pgx$Y)
      s1 <- grep("^[.]", px, value = TRUE, invert = TRUE)[1]
      shiny::updateSelectInput(session, "colorby", choices = px, selected = s1)
    })

    cor_scatter.DATA <- shiny::reactive({

      shiny::req(cor_gene, pgx$X, cor_table)
      
      this.gene <- cor_gene()
      NTOP <- 50
      R <- getGeneCorr()
      sel <- cor_table$rownames_current()
      sel <- head(intersect(sel, rownames(R)),NTOP)
      rho <- R[sel, "cor"]
      
      if (length(rho) == 1) names(rho) <- rownames(R)[1]
      pp <- unique(c(this.gene, names(rho)))
      X <- pgx$X[pp,]
      
      colorby <- input$colorby
      shiny::req(colorby)      
      pheno <- factor(pgx$samples[, colorby])
            
      dt <- list(
        rho = rho,
        X = X,
        this.gene = this.gene,
        pheno = pheno,
        colorby = colorby,
        COL = COL
      )
      
      return(dt)
    })

    plot_scatter <- function(nrow, ncol) {

      dt <- cor_scatter.DATA()
      shiny::req(dt)
      
      rho <- dt$rho
      X <- dt$X
      pheno <- dt$pheno
      colorby <- dt$colorby
      this.gene <- dt$this.gene
      COL <- rep(dt$COL,99)
      
      if (length(rho) == 0) {
        return(NULL)
      }
      
      klr <- COL[as.integer(pheno)]
      ndim <- length(pheno)

      cex_levels <- c(1.2, 0.8, 0.5, 0.2)*1.2
      dim_cuts <- c(0, 40, 100, 200, Inf)
      cex <- cex_levels[findInterval(ndim, dim_cuts)]

      nplots <- nrow*ncol
      rho <- head(rho, nplots)

      swapaxis <- input$swapaxis

      par(
        mfrow = c(nrow, ncol), mar = c(3, 3.5, 0.5, 0),
        mgp = c(1.7, 0.6, 0), oma = c(0, 0, 0.5, 0.5)
      )

      i <- 1
      for (i in 1:length(rho)) {
        gene2 <- names(rho)[i]
        if (swapaxis) {
          x <- X[gene2, ]
          y <- X[this.gene, ]
          xlab <- gene2
          ylab <- this.gene
        } else {
          y <- X[gene2, ]
          x <- X[this.gene, ]
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
          tt <- c("   ", levels(pheno))
          legend("topleft",
            legend = tt,
            fill = c(NA, COL), inset = c(0.02, 0.02),
            border = c(NA, rep("black", 99)),
            cex = 0.95, box.lwd = 0, pt.lwd = 0,
            x.intersp = 0.5, y.intersp = 0.8
          )
          legend("topleft", colorby,
            x.intersp = -0.2,
            cex = 0.95, y.intersp = 0.45, bty = "n"
          )
        }
      }
    }

    cor_scatter.PLOTFUN <- function() {
      if(input$layout=='3x3') {
        nrow = ncol = 3
      } else if(input$layout=='4x4') {
        nrow = ncol = 4
      } else {
        nrow = ncol = 5
      }
      plot_scatter(nrow, ncol)
    }

    cor_scatter.PLOTFUN2 <- function() {
      if(input$layout=='3x3') {
        nrow = 3
        ncol = 5
      } else if(input$layout=='4x4') {
        nrow = 4
        ncol = 6
      } else {
        nrow = 5
        ncol = 7
      }
      plot_scatter(nrow, ncol)
    }

    PlotModuleServer(
      "plot",
      plotlib = "base",
      func = cor_scatter.PLOTFUN,
      func2 = cor_scatter.PLOTFUN2,      
      csvFunc = cor_scatter.DATA, ##  NOTE: Not sure what should be the plot data!
      res = c(100, 120), ## resolution of plots
      pdf.width = 6,
      pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
