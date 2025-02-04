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
    info.methods,
    info.extra_link,
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
        selected = "3x3", inline = TRUE
      ),
      "Choose the layout for correlation plots.",
    ),
    withTooltip(
      shiny::checkboxInput(ns("swapaxis"), "Swap XY-axes"),
      "Transpose plot, i.e. swap X-axis and Y-axis.",
    )
  )

  PlotModuleUI(ns("plot"),
    title = title,
    plotlib = "plotly",
    label = "c",
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    caption = caption,
    options = cor_scatter.opts,
    download.fmt = c("png", "pdf"),
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
                                                sel_gene,
                                                COL,
                                                watermark = FALSE,
                                                labeltype) {
  moduleServer(id, function(input, output, session) {
    shiny::observe({
      px <- colnames(pgx$Y)
      s1 <- grep("^[.]", px, value = TRUE, invert = TRUE)[1]
      shiny::updateSelectInput(session, "colorby", choices = px, selected = s1)
    })

    cor_scatter.DATA <- shiny::reactive({
      shiny::req(sel_gene, pgx$X, cor_table)

      this.gene <- sel_gene()
      NTOP <- 50
      R <- getGeneCorr()
      sel <- cor_table$rownames_current()
      sel <- head(intersect(sel, rownames(R)), NTOP)
      shiny::req(sel)
      rho <- R[sel, "cor"]

      if (length(rho) == 1) names(rho) <- rownames(R)[1]
      pp <- unique(c(this.gene, names(rho)))
      X <- pgx$X[pp, , drop = FALSE]

      names(rho) <- playbase::probe2symbol(names(rho), pgx$genes, labeltype(), fill_na = TRUE)
      this.gene <- playbase::probe2symbol(this.gene, pgx$genes, labeltype(), fill_na = TRUE)
      rownames(X) <- playbase::probe2symbol(rownames(X), pgx$genes, labeltype(), fill_na = TRUE)

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

    plotly_scatter <- function(n_row, n_cols, markersize = 10, axis_title_pos = c(-0.07, -0.06),
                               margin_l = 50, margin_b = 10, interplot_margin = 0.03) {
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
      if (length(rho) < nplots) {
        n_row <- ceiling(sqrt(length(rho)))
      }

      # Assemble subplots
      sub_plots <- vector("list", length(rho))
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
          xlab <- this.gene
          ylab <- gene2
        }
        dx <- diff(range(x))
        dy <- diff(range(y))
        title_i <- gene2
        title_loc <- max(y) - 0.1 * dy
        # Make regression line
        fit <- lm(y ~ x)
        newdata <- data.frame(x = range(x))
        newdata$y <- predict(fit, newdata)
        plt <- plotly::plot_ly(
          hovertemplate = paste0(
            "<b>%{fullData.name}<br>",
            xlab,
            " Expression:</b> %{x}<br><b>",
            gene2,
            " Expression:</b> %{y}<extra></extra>"
          )
        ) %>%
          # Add the points
          plotly::add_trace(
            x = x,
            y = y,
            name = pheno,
            color = klr,
            type = "scatter",
            mode = "markers",
            marker = list(size = markersize),
            showlegend = (i == 1)
          ) %>%
          # Add the regression line
          plotly::add_trace(
            data = newdata,
            x = ~x,
            y = ~y,
            showlegend = FALSE,
            type = "scatter",
            line = list(color = "rgb(22, 96, 167)", dash = "dot"),
            mode = "lines",
            inherit = FALSE
          ) %>%
          # Legend
          plotly::layout(
            legend = list(orientation = "h", bgcolor = "transparent"),
            xaxis = list(
              title = list(text = xlab, font = list(size = 11), standoff = 5),
              range = list(min(x) - 0.05 * dx, max(x) + 0.05 * dx)
            ),
            yaxis = list(
              title = list(text = ylab, font = list(size = 11), standoff = 0),
              range = list(min(y) - 0.05 * dx, max(y) + 0.05 * dy)
            )
          ) %>%
          # Plot title
          plotly::add_annotations(
            # text = paste("<b>", title_i, "</b>"),
            text = "",
            font = list(size = 10),
            showarrow = FALSE,
            xanchor = "left",
            yanchor = "bottom",
            x = min(x),
            y = title_loc
          )
        if (pgx$datatype != "scRNAseq") {
          plt %>% playbase::plotly_build_light(.)
        }
        sub_plots[[i]] <- plt
      }

      # Assemble all subplot in to grid
      suppressWarnings(
        all_plt <- plotly::subplot(
          sub_plots,
          nrows = n_row,
          margin = interplot_margin,
          titleY = TRUE,
          titleX = TRUE
        ) %>%
          # Add common axis titles
          plotly::layout(
            annotations = list(
              list(
                x = 0.5,
                y = axis_title_pos[2],
                text = glue::glue(tspan("gene expression (log2)", js = FALSE)),
                # text = xlab,
                font = list(size = 15),
                showarrow = FALSE,
                xref = "paper",
                yref = "paper"
              ),
              list(
                y = 0.5,
                x = axis_title_pos[1],
                text = glue::glue(tspan("gene expression (log2)", js = FALSE)),
                # text = xlab,
                textangle = 270,
                font = list(size = 15),
                showarrow = FALSE,
                xref = "paper",
                yref = "paper"
              )
            ),
            margin = list(
              l = margin_l,
              b = margin_b
            )
          )
      )
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
      plotly_scatter(
        nrow, ncol,
        markersize = 6,
        axis_title_pos = c(-0.12, -0.1),
        interplot_margin = 0.04,
        margin_l = 65,
        margin_b = 60
      )
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
      plotly_scatter(
        nrow, ncol,
        markersize = 8,
        axis_title_pos = c(-0.05, -0.11),
        interplot_margin = 0.045,
        margin_l = 90,
        margin_b = 80
      )
    }

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = cor_scatter.PLOTFUN,
      func2 = cor_scatter.PLOTFUN2,
      res = c(100, 120), ## resolution of plots
      pdf.width = 6,
      pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
