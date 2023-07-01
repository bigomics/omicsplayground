##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

dataview_plot_abundance_ui <- function(
  id,
  label = "",
  height,
  width,
  title,
  info.text,
  caption
  ) {
  ns <- shiny::NS(id)

  menu_grouped <- "<code>grouped</code>"
  info_text <- paste0("Barplot showing the percentage of counts in terms of major gene types such as ribosomal protein genes, kinases or RNA binding motifs for each group. The samples (or cells) can be grouped/ungrouped in the ", menu_grouped, " setting uder the main <i>Options</i>.")

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    caption = caption,
    plotlib = "plotly",
    info.text = info_text,
    options = NULL,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height
  )
}

dataview_plot_abundance_server <- function(id,
                                           getCountsTable,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## extract data from pgx object
    plot_data <- shiny::reactive({
      res <- getCountsTable()
      shiny::req(res)
      res <- list(
        prop.counts = res$prop.counts
      )
      res
    })

    ## plot.RENDER <- function() {



    ##   klr <- colorRampPalette(
    ##     c(


    ##     ),

    ##   )(nrow(res$prop.counts))



    ##   if (length(names.arg) > 20) {

    ##   }



    ##   barplot(res$prop.counts,

    ##     cex.lab = 1.0, border = NA,
    ##     ylim = c(0, ymax) * 1.6, ylab = "abundance (%)",
    ##     names.arg = names.arg, cex.names = cex.names,

    ##   )
    ##   leg <- legend("topleft",

    ##     fill = rev(klr), cex = 1, y.intersp = 0.75, bty = "n", plot = FALSE
    ##   )




    ##   legend(
    ##     x = c(leftx, rightx), y = c(topy, bottomy),

    ##     fill = rev(klr), bty = "n", cex = 0.9, y.intersp = 0.75
    ##   )
    ## }
    ## 
    ## modal_plot.RENDER <- function() {

    ## }

    plotly.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)

      long.data <- reshape2::melt(head(res$prop.counts, 5))
      colnames(long.data) <- c("gene", "sample", "value")

      ## stacked barchart
      fig <-
        plotly::plot_ly(
          data = long.data,
          x = ~sample,
          y = ~value,
          type = "bar",
          color = ~gene,
          colors = omics_pal_d("muted")(length(unique(long.data$gene))),
          hovertemplate = ~ paste0(
            "Sample: <b>", sample, "</b><br>",
            "Gene: <b>", gene, "</b><br>",
            "Cum. proportion: <b>", sprintf("%2.1f", value), "%</b>",
            "<extra></extra>"
          )
        ) %>%
        plotly::layout(
          barmode = "stack",
          xaxis = list(title = FALSE),
          yaxis = list(title = "Cumulative proportion", ticksuffix = "%"),
          font = list(family = "Lato"),
          margin = list(l = 10, r = 10, b = 10, t = 10)
        ) %>%
        plotly_default()
      fig
    }

    modal_plotly.RENDER <- function() {
      fig <- plotly.RENDER() %>%
        plotly_modal_default() %>%
        plotly::layout(
          showlegend = FALSE
        )  
      fig
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",

      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(90, 170), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
