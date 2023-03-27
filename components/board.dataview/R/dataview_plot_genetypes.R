##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

dataview_plot_genetypes_ui <- function(id, label = "", height = c(600, 800)) {
  ns <- shiny::NS(id)

  menu_grouped <- "<code>grouped</code>"
  info_text <- paste0("Barplot showing the percentage of counts in terms of major gene types such as ribosomal protein genes, kinases or RNA binding motifs for each group. The samples (or cells) can be grouped/ungrouped in the ", menu_grouped, " setting uder the main <i>Options</i>.")

  PlotModuleUI(
    ns("pltmod"),
    title = "Relative abundance of major gene types",
    label = label,
    plotlib = "plotly",
    info.text = info_text,
    options = NULL,
    download.fmt = c("png", "pdf", "csv"),
    width = c("auto", "100%"),
    height = height
  )
}

dataview_plot_genetypes_server <- function(id,
                                           getCountsTable,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## extract data from pgx object
    plot_data <- shiny::reactive({
      res <- getCountsTable()
      shiny::req(res)
      res <- list(
        prop.counts = res$prop.counts,
        gset.genes = res$gset.genes
      )
      res
    })

    plot.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)

      nsamples <- ncol(res$prop.counts)
      klr <- colorRampPalette(
        c(
          rgb(0.2, 0.5, 0.8, 0.8),
          rgb(0.2, 0.5, 0.8, 0.1)
        ),
        alpha = TRUE
      )(nsamples)

      ymax <- max(colSums(res$prop.counts, na.rm = TRUE))
      names.arg <- colnames(res$prop.counts)
      if (length(names.arg) > 20) {
        names.arg <- rep("", length(names.arg))
      }
      cex.names <- ifelse(length(names.arg) > 10, 0.8, 0.9)

      par(mfrow = c(1, 2), mar = c(4, 0, 2, 0.5), mgp = c(2.2, 0.8, 0))
      frame()
      barplot(
        t(res$prop.counts) / nsamples,
        horiz = TRUE, las = 1,
        cex.lab = 1.0, border = NA,
        # ylim = c(0,ymax)*1.6, ylab = "abundance (%)",
        ## ylim = c(0,ymax)*1.6,
        xlab = "abundance (%)",
        # names.arg = names.arg, cex.names = cex.names,
        col = klr
      )

      leg <- legend("topright",
        legend = rev(colnames(res$prop.counts)),
        fill = rev(klr), cex = 1, y.intersp = 0.75, bty = "n", plot = FALSE
      )
      leftx <- leg$rect$left * 0.9
      rightx <- leg$rect$right * 0.9
      topy <- leg$rect$top
      bottomy <- leg$rect$bottom
      legend(
        x = c(leftx, rightx), y = c(topy, bottomy),
        legend = rev(colnames(res$prop.counts)),
        fill = rev(klr), bty = "n", cex = 0.9, y.intersp = 0.75
      )
    }

    modal_plot.RENDER <- function() {
      plot.RENDER()
    }

    plotly.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)

      avg.prop <- head(rowMeans(res$prop.counts), 15)
      genes <- head(res$gset.genes, 15)
      family <- paste0(names(avg.prop), "  ")
      family <- factor(family, levels = family)

      df <- data.frame(family = family, prop = avg.prop, genes = genes)

      ## stacked barchart
      fig <-
        plotly::plot_ly(
          data = df,
          x = ~prop,
          y = ~family,
          type = "bar",
          marker = list(
            color = omics_colors("brand_blue")
          ),
          hovertemplate = ~ paste0(
            "Gene family: <b>", family, "</b><br>",
            ## NOTE: tooltip looks awful due to way too many genes
            ## TODO: discuss potential solutions; how about showing number of genes or the top 3?
            "Genes: <b>", genes, "</b>",
            "<extra></extra>"
          )
        ) %>%
        plotly::layout(
          yaxis = list(title = FALSE),
          xaxis = list(title = "Proportion", ticksuffix = "%"),
          font = list(family = "Lato"),
          margin = list(l = 10, r = 10, b = 10, t = 10),
          showlegend = FALSE
        ) %>%
        plotly_default()
      fig
    }

    modal_plotly.RENDER <- function() {
      fig <- plotly.RENDER() %>%
        plotly_modal_default() %>%
        plotly::layout(
          showlegend = TRUE
        )
      ## fig <- plotly::style(fig, marker.size = 14)
      fig
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      ## plotlib2 = "base",
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(90, 170), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
