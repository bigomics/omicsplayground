##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' Word Cloud Plot UI
#'
#' @description
#' Creates the UI for the word cloud plot module.
#'
#' @param id Module ID string
#' @param title Plot title
#' @param label Plot label
#' @param info.text Info text
#' @param caption Caption text
#' @param height Plot height
#' @param width Plot width
#'
#' @return
#' A Shiny Module UI definition
wordcloud_plot_wordcloud_ui <- function(
    id,
    title,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    caption,
    height) {
  ns <- shiny::NS(id)

  info_text <- "<strong>Word cloud.</strong> Word cloud of the most enriched keywords for the data set. Select a keyword in the 'Enrichment table'. In the plot settings, users can exclude certain words from the figure, or choose the color palette. The sizes of the words are relative to the normalized enrichment score (NES) from the GSEA computation. Keyword enrichment is computed by running GSEA on the mean (squared) enrichment profile (averaged over all contrasts). For each keyword, we defined the 'keyword set' as the collection of genesets that contain that keyword in the title/description."

  wordcloud_opts <- shiny::tagList(
    withTooltip(shiny::selectInput(ns("wordcloud_exclude"), "Exclude words:", choices = NULL, multiple = TRUE),
      "Paste a keyword to exclude it from the plot.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(
      shiny::selectInput(ns("wordcloud_colors"), "Colors:",
        choices = c("Blues", "Greys", "Accent", "Dark2"),
        multiple = FALSE
      ),
      "Choose a set of colors.",
      placement = "top", options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = "b",
    info.text = info_text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    options = wordcloud_opts,
    height = height,
    download.fmt = c("png", "pdf", "csv", "svg")
  )
}


#' Word Cloud Plot Server Function
#'
#' @description
#' Server function for rendering the word cloud plot in the wordcloud module.
#'
#' @param id Module id string.
#' @param getCurrentWordEnrichment Reactive returning the current word enrichment data.
#' @param watermark Watermark text.
#'
#' @details
#' This function generates the word cloud plot based on the reactive
#' getCurrentWordEnrichment. It handles rendering the plot, updating based on
#' input settings, and adding the watermark.
#'
#' @return
#' Shiny module server function.
wordcloud_plot_wordcloud_server <- function(id,
                                            getCurrentWordEnrichment,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      getCurrentWordEnrichment()
    })

    wordcloud.RENDER <- shiny::reactive({
      df <- plot_data()

      ## update selectors
      words <- sort(df$word)
      shiny::updateSelectInput(session, "wordcloud_exclude",
        choices = words,
        selected = input$wordcloud_exclude
      )

      excl.words <- input$wordcloud_exclude

      if (length(excl.words) > 0) {
        df <- df[which(!df$word %in% excl.words), ]
      }

      cex1 <- 1 + round((5 * rank(abs(df$NES)) / nrow(df))**2)
      cex2 <- (-log10(df$padj))**1.0
      size <- 10 * abs(cex1 * cex2)**1
      size[is.na(size)] <- 0
      minsize <- tail(sort(size), 250)[1]

      color.pal <- input$wordcloud_colors

      par(mar = c(1, 1, 1, 1) * 0)
      suppressWarnings(suppressMessages(
        wordcloud::wordcloud(
          words = df$word, freq = size,
          colors = RColorBrewer::brewer.pal(8, color.pal),
          scale = c(2, 0.1) * 0.9, min.freq = minsize
        )
      ))
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = wordcloud.RENDER,
      csvFunc = plot_data,
      pdf.width = 5, pdf.height = 5,
      res = 72,
      add.watermark = watermark
    )
  })
}
