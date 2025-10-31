##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wordcloud_plot_wordtsne_ui <- function(
  id,
  height,
  title,
  info.text,
  info.methods,
  info.references,
  info.extra_link,
  caption
) {
  ns <- shiny::NS(id)

  wordtsne_options <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(ns("tsne_algo"), "Clustering algorithm:",
        choices = c("tsne", "umap"), inline = TRUE
      ),
      "Choose a clustering algorithm: t-SNE or UMAP."
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = "c",
    plotlib = "plotly",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    options = wordtsne_options,
    height = height,
    download.fmt = c("png", "pdf", "svg")
  )
}

wordcloud_plot_wordtsne_server <- function(id,
                                           getCurrentWordEnrichment,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    wordtsne.PLOTLY <- shiny::reactive({
      topFreq <- getCurrentWordEnrichment()

      df <- topFreq
      klr <- ifelse(df$padj <= 0.05, "red", "grey")
      ps1 <- 0.5 + 3 * (1 - df$padj) * (df$NES / max(df$NES))**3

      ## label top 20 words
      df$label <- rep("", nrow(df))
      jj <- head(order(-abs(df$NES)), 20)
      df$label[jj] <- as.character(df$word[jj])
      cex <- 1
      #
      df$abs.NES <- abs(df$NES)**2

      if (input$tsne_algo == "tsne") {
        pos <- cbind(x = df$tsne.x, y = df$tsne.y)
      } else {
        pos <- cbind(x = df$umap.x, y = df$umap.y)
      }

      plt <- plotly::plot_ly(
        df,
        text = df$word, hoverinfo = "text"
        #
      ) %>%
        plotly::add_markers(
          type = "scatter",
          x = pos[, 1], y = pos[, 2],
          color = klr,
          size = ~abs.NES,
          marker = list(
            line = list(color = "grey20", width = 0.6)
          )
        ) %>%
        plotly::add_annotations(
          x = pos[, 1], y = pos[, 2],
          text = df$label,
          font = list(size = 12),
          showarrow = FALSE
        )

      ax <- list(
        title = "",
        showticklabels = FALSE,
        showgrid = FALSE
      )

      m <- list(
        l = 0,
        r = 0,
        b = 0,
        t = 0,
        pad = 4
      )

      plt <- plt %>%
        plotly::layout(
          xaxis = ax,
          yaxis = ax,
          showlegend = FALSE,
          margin = m
        )
      return(plt)
    })

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = wordtsne.PLOTLY,
      pdf.width = 5, pdf.height = 5,
      res = 72,
      add.watermark = watermark
    )
  })
}
