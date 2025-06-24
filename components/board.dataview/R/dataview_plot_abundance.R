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
    caption) {
  ns <- shiny::NS(id)

  menu_grouped <- "<code>grouped</code>"

  options <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("show_overall_prop"), "Show overall proportions", FALSE
      ),
      "Show proportions of major feature (gene, protein) types in the whole dataset.",
      placement = "top"
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    caption = caption,
    plotlib = "plotly",
    info.text = info.text,
    options = options,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height
  )
}

dataview_plot_abundance_server <- function(id,
                                           getCountsTable,
                                           r.samples = reactive(""),
                                           watermark = FALSE) {

  moduleServer(id, function(input, output, session) {

    plot_data <- shiny::reactive({
      res <- getCountsTable()
      samples <- r.samples()
      shiny::req(res)
      if (!input$show_overall_prop) {
        return(list(prop.counts = res$prop.counts))
      } else {
        return(list(prop.counts = res$prop.counts, gset.genes = res$gset.genes))
      }
    })

    plotly.RENDER <- function(return_csv = FALSE) {
      res <- plot_data()
      shiny::req(res)
      
      if (!input$show_overall_prop) {

        long.data <- reshape2::melt(head(res$prop.counts, 5))
        colnames(long.data) <- c("gene", "sample", "value")
        if (return_csv) return(long.data)
        
        ## stacked barchart
        fig <-
          plotly::plot_ly(
            data = long.data,
            x = ~sample,
            y = ~value,
            type = "bar",
            color = ~gene,
            colors = omics_pal_d(palette = "expanded")(length(unique(long.data$gene))),
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
      } else {

        avg.prop <- head(rowMeans(res$prop.counts, na.rm = TRUE), 15)
        genes <- head(res$gset.genes, 15)
        family <- paste0(names(avg.prop), "  ")
        family <- factor(family, levels = family)
        df <- data.frame(family = family, prop = avg.prop, genes = genes)
        if (return_csv) return(df)

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

    }
    
    modal_plotly.RENDER <- function() {
      fig <- plotly.RENDER() %>%
        plotly_modal_default() %>%
        plotly::layout(showlegend = FALSE)
      fig
    }

    plot_data_csv <- function() {
      df <- plotly.RENDER(return_csv = TRUE)
      return(df)
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      csvFunc = plot_data_csv,
      res = c(90, 170),
      pdf.width = 6,
      pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
