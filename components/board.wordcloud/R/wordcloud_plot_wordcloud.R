##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

wordcloud_plot_wordcloud_ui <- function(id, height) {
  ns <- shiny::NS(id)

  info_text <- "<strong>Word cloud.</strong> Word cloud of the most enriched keywords for the data set. Select a keyword in the 'Enrichment table'. In the plot settings, users can exclude certain words from the figure, or choose the color palette. The sizes of the words are relative to the normalized enrichment score (NES) from the GSEA computation. Keyword enrichment is computed by running GSEA on the mean (squared) enrichment profile (averaged over all contrasts). For each keyword, we defined the 'keyword set' as the collection of genesets that contain that keyword in the title/description."

  wordcloud_opts = shiny::tagList(
    withTooltip(shiny::selectInput(ns("wordcloud_exclude"),"Exclude words:", choices=NULL, multiple=TRUE),
                "Paste a keyword to exclude it from the plot.", placement="top", options = list(container = "body")),
    withTooltip(shiny::selectInput(ns("wordcloud_colors"),"Colors:", choices=c("Blues","Greys","Accent","Dark2"),
                                   multiple=FALSE),
                "Choose a set of colors.", placement="top", options = list(container = "body"))
  )

  PlotModuleUI(
    ns("plot"),
    title = "Word cloud",
    label = "b",
    info.text = info_text,
    options = wordcloud_opts,
    height = height,
    download.fmt = c("png", "pdf")
  )
}

wordcloud_plot_wordcloud_server <- function(id,
                                            getCurrentWordEnrichment,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    wordcloud.RENDER <- shiny::reactive({

      topFreq <- getCurrentWordEnrichment()
      df <- topFreq

      ## update selectors
      words <- sort(df$word)
      shiny::updateSelectInput(session, "wordcloud_exclude", choices=words,
                               selected = input$wordcloud_exclude)

      excl.words <- input$wordcloud_exclude

      if(length(excl.words)>0) {
        df <- df[ which(!df$word %in% excl.words), ]
      }

      cex1 <- 1+round((5*rank(abs(df$NES))/nrow(df))**2)
      cex2 <- (-log10(df$padj))**1.0
      size <- 10*abs(cex1 * cex2)**1
      minsize <- tail(sort(size),250)[1]

      color.pal = input$wordcloud_colors

      par(mar=c(1,1,1,1)*0)
      suppressWarnings(suppressMessages(
        wordcloud::wordcloud(
          words = df$word, freq = size,
          colors = RColorBrewer::brewer.pal(8, color.pal),
          scale=c(2,0.1)*0.9, min.freq=minsize)
      ))
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = wordcloud.RENDER,
      pdf.width = 5, pdf.height = 5,
      res=72,
      add.watermark = watermark
    )
  })
}
