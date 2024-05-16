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
        ns("colorby"), "Color by:",
        choices = NULL, multiple = FALSE
      ),
      "Variable to split and color by groups.",
      placement = "top"
    ),
    withTooltip(
      shiny::radioButtons(ns("layout"), "Layout:", c("3x3", "4x4", "5x5"),
        selected = "4x4", inline = TRUE
      ),
      "Choose the layout for correlation plots.",
    ),
    withTooltip(
      shiny::checkboxInput(ns("swapaxis"), "Swap XY-axes"),
      "Transpose plot, i.e. swap X-axis and Y-axis.",
    )
  )

  PlotModuleUI(ns("pltmod"),
    title = title,
    plotlib = "plotly",
    label = "c",
    info.text = info.text,
    caption = caption,
    options = cor_scatter.opts,
    download.fmt = c("png", "pdf"),
    width = width,
    height = height,
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
      shiny::updateSelectInput(session, "pltmod-subplot_selector", choices = subplot_names())
    })

    cor_scatter.DATA <- shiny::reactive({
      shiny::req(cor_gene, pgx$X, cor_table)

      this.gene <- cor_gene()
      NTOP <- 50
      R <- getGeneCorr()
      sel <- cor_table$rownames_current()
      sel <- head(intersect(sel, rownames(R)), NTOP)
      rho <- R[sel, "cor"]

      if (length(rho) == 1) names(rho) <- rownames(R)[1]
      pp <- unique(c(this.gene, names(rho)))
      X <- pgx$X[pp, ]

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

    subplot_names <- shiny::reactive({
      shiny::req(cor_gene(), pgx$X, cor_scatter.DATA(), input$layout)
      dt <- cor_scatter.DATA()
      rho <- dt$rho
      shiny::req(length(rho) > 0)
      if (input$layout == "3x3") {
        nrow <- ncol <- 3
      } else if (input$layout == "4x4") {
        nrow <- ncol <- 4
      } else {
        nrow <- ncol <- 5
      }
      nplots <- nrow * ncol
      names_download <- names(rho)[1:nplots]
      download_options <- 1:nplots
      names(download_options) <- names_download
      download_options <- c("All", download_options)
      return(download_options)
    })

    plotly_scatter <- function(n_row, n_cols, markersize = 10, axis_title_pos = c(-0.1, -0.1),
                               margin_l = 50, margin_b = 10, interplot_margin = 0.02) {
      # Load input data
      dt <- cor_scatter.DATA()

      shiny::req(dt)
      swapaxis <- input$swapaxis
      rho <- dt$rho
      X <- dt$X
      pheno <- dt$pheno
      colorby <- dt$colorby
      this.gene <- dt$this.gene
      COL <- rep(dt$COL, 99)

      shiny::req(length(rho) > 0)
      klr <- COL[as.integer(pheno)]

      nplots <- n_row * n_cols
      rho <- head(rho, nplots)
      # Aseemble subplots
      sub_plots <- vector("list", length(rho))
      for (i in 1:nplots) {
        gene2 <- names(rho)[i]
        if (swapaxis) {
          x <- X[gene2, ]
          y <- X[this.gene, ]
          xlab <- "Y"
          ylab <- this.gene
        } else {
          y <- X[gene2, ]
          x <- X[this.gene, ]
          xlab <- this.gene
          ylab <- "Y"
        }
        title_i <- gene2
        title_loc <- max(y) + 0.05 * max(y)
        # Make regression line
        fit <- lm(y ~ x)
        newdata <- data.frame(x = range(x))
        newdata$y <- predict(fit, newdata)
        plt <- plotly::plot_ly() %>%
          # Add the points
          plotly::add_trace(
            x = x, y = y, name = pheno, color = klr, type = "scatter", mode = "markers",
            marker = list(size = markersize),
            showlegend = i == 14
          ) %>%
          # Add the regression line
          plotly::add_trace(
            data = newdata,
            x = ~x, y = ~y, showlegend = FALSE, type = "scatter",
            line = list(color = "rgb(22, 96, 167)", dash = "dot"), mode = "lines",
            inherit = FALSE
          ) %>%
          # Legend
          plotly::layout(
            legend = list(orientation = "h", bgcolor = "transparent")
          ) %>%
          # Plot title
          plotly::add_annotations(
            text = paste("<b>", title_i, "</b>"),
            font = list(size = 10),
            showarrow = FALSE,
            xanchor = "left",
            yanchor = "bottom",
            x = 0,
            y = title_loc
          ) %>% playbase::plotly_build_light(.)
        sub_plots[[i]] <- plt
      }


      # Assemble all subplot in to grid
      suppressWarnings(
        all_plt <- plotly::subplot(sub_plots,
          nrows = n_row, margin = interplot_margin,
          titleY = FALSE, titleX = FALSE
        ) %>%
          # Add common axis titles
          plotly::layout(
            annotations = list(
              list(
                x = axis_title_pos[1], y = 0.5, text = glue::glue("<b> Gene {ylab} Expression </b>"),
                font = list(size = 13),
                textangle = 270,
                showarrow = FALSE, xref = "paper", yref = "paper"
              ),
              list(
                x = 0.5, y = axis_title_pos[2], text = glue::glue("<b> Gene {xlab} expression </b>"),
                font = list(size = 13),
                showarrow = FALSE, xref = "paper", yref = "paper"
              )
            ),
            margin = list(l = margin_l, b = margin_b)
          )
      )
      all_plt$plots <- sub_plots
      return(all_plt)
    }


    cor_scatter.PLOTFUN <- function() {
      if (input$layout == "3x3") {
        nrow <- ncol <- 3
      } else if (input$layout == "4x4") {
        nrow <- ncol <- 4
      } else {
        nrow <- ncol <- 5
      }
      plotly_scatter(nrow, ncol, markersize = 5, margin_l = 50)
    }

    cor_scatter.PLOTFUN2 <- function() {
      if (input$layout == "3x3") {
        nrow <- 3
        ncol <- 5
      } else if (input$layout == "4x4") {
        nrow <- 4
        ncol <- 6
      } else {
        nrow <- 5
        ncol <- 7
      }
      plotly_scatter(nrow, ncol, markersize = 10, axis_title_pos = c(-0.05, -0.1))
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = cor_scatter.PLOTFUN,
      func2 = cor_scatter.PLOTFUN2,
      res = c(100, 120), ## resolution of plots
      pdf.width = 6,
      pdf.height = 6,
      add.watermark = watermark,
      subplot = TRUE
    )
  }) ## end of moduleServer
}
