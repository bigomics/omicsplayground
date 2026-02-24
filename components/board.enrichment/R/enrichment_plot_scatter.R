##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' Enrichment Scatter Plot UI
#'
#' @description
#' Creates the UI for the enrichment scatter plot module.
#'
#' @param id Module ID string
#' @param title Plot title
#' @param label Plot label
#' @param info.text Info text to be displayed
#' @param caption Caption text
#' @param height Plot height
#' @param width Plot width
#'
#' @return
#' A Shiny Module UI definition
enrichment_plot_scatter_ui <- function(
  id,
  title,
  info.text,
  info.methods,
  info.extra_link,
  caption,
  height,
  width
) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    plotlib = "plotly",
    title = title,
    label = "d",
    caption = caption,
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg"),
    editor = TRUE,
    ns_parent = ns,
    plot_type = "clustering",
    palette_default = "original"
  )
}


#' Enrichment Scatter Plot Server Function
#'
#' @description Server function for generating an enrichment analysis
#' scatter plot in a Shiny app.
#'
#' @param id Shiny module id
#' @param pgx PGX object
#' @param getGSEAReactive Reactive expression for getting GSEA results
#' @param getGeneSelected Reactive expression for selected gene
#' @param getGsetSelected Reactive expression for selected gene set
#' @param watermark Add watermark to plot? Logical
#'
#' @return None. Generates scatter plot.
enrichment_plot_scatter_server <- function(id,
                                           pgx,
                                           gene_selected,
                                           gs_contrast,
                                           subplot.MAR,
                                           gset_selected,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    getcolors <- function(pgx, comp0) {
      ## get colors (what a mess...)
      contr.matrix <- pgx$model.parameters$contr.matrix

      exp.matrix <- pgx$model.parameters$exp.matrix
      xgroup <- as.character(pgx$Y$group)

      grp.name <- strsplit(comp0, split = "[._ ]vs[._ ]")[[1]]
      grp.name <- c(grp.name, "other")
      xsign <- sign(exp.matrix[, comp0])
      xgroup <- grp.name[1 * (xsign > 0) + 2 * (xsign < 0) + 1 * (xsign == 0)]

      names(xgroup) <- rownames(pgx$Y)
      samples <- names(which(exp.matrix[, comp0] != 0))

      xgroup1 <- xgroup[samples]
      ngrp <- length(unique(xgroup1))
      grp.klr <- c("grey90", rep(RColorBrewer::brewer.pal(12, "Paired"), 99)[1:ngrp])
      names(grp.klr) <- c("other", as.character(sort(unique(xgroup1))))

      xgroup2 <- xgroup
      xgroup2[which(!(xgroup %in% xgroup1))] <- "other"
      sample.klr <- grp.klr[as.character(xgroup2)]
      names(sample.klr) <- rownames(pgx$samples)
      list(samples = sample.klr, group = grp.klr, xgroup = xgroup2)
    }

    basesubplot_scatter.RENDER <- function() {
      par(mfrow = c(1, 1), mgp = c(1.8, 0.8, 0), oma = c(0, 0, 0, 0.4))
      par(mar = subplot.MAR)
      shiny::req(pgx$X)
      gene <- rownames(pgx$X)[1]
      sel <- gene_selected()
      gset <- gset_selected()
      shiny::req(sel, gset)

      not.selected <- (is.null(sel) || length(sel) == 0)
      shiny::validate(shiny::need(
        not.selected == FALSE, tspan("Please select a gene", js = FALSE)
      ))

      gene <- sel$gene
      gset <- gset[1]
      gx <- pgx$X[sel$rn, ]
      sx <- pgx$gsetX[gset, ]
      if (length(gx) == 0 || length(sx) == 0 ||
        length(gx) != length(sx)) {
        frame()
        return(NULL)
      }
      ## get colors
      comp0 <- gs_contrast()
      klrs <- getcolors(pgx, comp0)
      klr <- klrs$samples[names(sx)]
      klr <- paste0(gplots::col2hex(klr), "99")

      cex1 <- c(1.4, 0.8, 0.3)[cut(length(gx), c(0, 100, 500, 99999))]
      gset1 <- playbase::breakstring(substring(gset, 1, 80), 32)
      tt <- paste(playbase::breakstring(gset, 40, 80), " vs. ", gene)
      base::plot(gx, sx,
        col = klr, main = tt,
        ylab = tspan("gene set enrichment", js = FALSE),
        xlab = paste(gene, "expression"),
        cex.lab = 1, pch = 19, cex = 1.0 * cex1, cex.main = 0.85
      )
      abline(lm(sx ~ gx), lty = 2, lwd = 0.7, col = "black")
    }

    ## Editor: dynamic color pickers for custom palette
    output$custom_palette_ui <- shiny::renderUI({
      shiny::req(input$palette == "custom")
      comp0 <- gs_contrast()
      shiny::req(comp0)
      klrs <- getcolors(pgx, comp0)
      groups <- names(klrs$group)
      shiny::req(length(groups) > 0)
      pickers <- lapply(seq_along(groups), function(i) {
        colourpicker::colourInput(
          session$ns(paste0("custom_color_", i)),
          label = groups[i],
          value = klrs$group[i]
        )
      })
      shiny::tagList(pickers)
    })

    plotlysubplot_scatter <- function() {
      shiny::req(pgx$X)
      gene <- rownames(pgx$X)[1]
      sel <- gene_selected()
      gset <- gset_selected()
      shiny::req(sel, gset)
      gene <- sel$gene
      gene <- playbase::probe2symbol(gene, pgx$genes, "gene_name", fill_na = TRUE)
      gset <- gset[1]
      gx <- pgx$X[sel$rn, ]
      sx <- pgx$gsetX[gset, ]
      if (length(gx) == 0 || length(sx) == 0 ||
        length(gx) != length(sx)) {
        frame()
        return(NULL)
      }
      ## get colors and group assignments
      comp0 <- gs_contrast()
      klrs <- getcolors(pgx, comp0)
      xgroup <- klrs$xgroup[names(sx)]
      grp.klr <- klrs$group
      groups <- names(grp.klr)
      clrs.length <- length(groups)

      ## Editor: palette override
      palette <- input$palette
      if (!is.null(palette) && palette == "custom") {
        COL <- sapply(seq_len(clrs.length), function(j) {
          val <- input[[paste0("custom_color_", j)]]
          if (is.null(val)) grp.klr[j] else val
        })
        names(COL) <- groups
      } else if (!is.null(palette) && palette != "original") {
        pal_colors <- rep(omics_pal_d(palette = palette)(8), ceiling(clrs.length / 8))
        COL <- c("grey90", pal_colors[1:(clrs.length - 1)])
        names(COL) <- groups
      } else {
        COL <- grp.klr
      }

      pheno <- factor(xgroup, levels = groups)
      gset1 <- playbase::breakstring(substring(gset, 1, 80), 32)
      tt <- paste(playbase::breakstring(gset, 40, 80), " vs. ", gene)

      # Assemble Plot
      fit <- lm(sx ~ gx)
      newdata <- data.frame(gx = range(gx))
      newdata$sx <- predict(fit, newdata)
      plt <- plotly::plot_ly(
        hovertemplate = paste0(
          "<b>%{fullData.name}<br>",
          gene,
          " Expression:</b> %{x}<br><b>",
          tspan("Geneset", js = FALSE),
          " Enrichment:</b> %{y}<extra></extra>"
        )
      ) %>%
        # Axis
        plotly::layout(
          xaxis = list(title = paste(gene, "expression"), titlefont = 5),
          yaxis = list(title = tspan("gene set enrichment", js = FALSE), titlefont = 5),
          legend = list(
            x = 0.05, y = 1.1, xanchor = "center",
            orientation = "h", bgcolor = "transparent",
            font = list(size = 7)
          )
        ) %>%
        # Add the points
        plotly::add_trace(
          x = gx, y = sx, color = pheno, colors = COL,
          type = "scatter", mode = "markers",
          marker = list(size = 5),
          showlegend = TRUE
        ) %>%
        # Add the regression line
        plotly::add_trace(
          data = newdata,
          x = ~gx, y = ~sx, showlegend = FALSE, type = "scatter",
          line = list(color = "rgb(22, 96, 167)", dash = "dot"), mode = "lines",
          inherit = FALSE
        )

      return(plt)
    }

    subplotly_scatter.RENDER <- list(
      card = function() {
        plotlysubplot_scatter() %>% plotly_default()
      },
      expand = function() {
        plotlysubplot_scatter() %>% plotly_default()
      }
    )

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = subplotly_scatter.RENDER$card, # basesubplot_scatter.RENDER,
      func2 = subplotly_scatter.RENDER$expand,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 100),
      add.watermark = watermark,
      parent_session = session
    )
  })
}
