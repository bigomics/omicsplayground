##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

WordCloudInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    withTooltip(shiny::selectInput(ns("wc_contrast"), "Contrast:", choices = NULL),
      "Select the contrast corresponding to the comparison of interest.",
      placement = "top"
    )
  )
}

WordCloudUI <- function(id) {
  fullH <- 750
  rowH <- 660 ## row height of panel
  tabH <- 200 ## row height of panel
  tabH <- "70vh" ## row height of panel
  halfH <- c("calc(50vh - 70px)",TABLE_HEIGHT_MODAL)
  
  ns <- shiny::NS(id) ## namespace
  shiny::tabsetPanel(
    id = ns("tabs"),
    tabs <- shiny::tabPanel(
      "",
      div(
        class = "row",
        div(
          class = "col-md-4",
          wordcloud_plot_enrichment_ui(
            ns("gseaplots"), 
            title = "Enrichment plots",
            info.text = "Select a keyword by clicking a word in the 'Enrichment table'. Keyword enrichment is computed by running GSEA on the enrichment score profile for all contrasts. We defined the test set as the collection of genesets that contain the keyword in the title/description. Black vertical bars indicate the position of gene sets that contains the *keyword* in the ranked list of enrichment scores. The curve in green corresponds to the 'running statistic' of the keyword enrichment score. The more the green ES curve is shifted to the upper left of the graph, the more the keyword is enriched in the first group. Conversely, a shift of the green ES curve to the lower right, corresponds to keyword enrichment in the second group.",
            caption = "Keyword enrichment plots for all available contrasts. ",
            height = halfH
            )
        ),
        div(
          class = "col-md-4",
          wordcloud_plot_wordcloud_ui(
            ns("wordcloud"),
            height = halfH,
            title = "Word cloud",
            info.text = "In the plot settings, users can exclude certain words from the figure, or choose the color palette. The sizes of the words are relative to the normalized enrichment score (NES) from the GSEA computation. Keyword enrichment is computed by running GSEA on the mean (squared) enrichment profile (averaged over all contrasts). For each keyword, we defined the 'keyword set' as the collection of genesets that contain that keyword in the title/description.",
            caption = "Word cloud of the most enriched keywords for the data set. Words are taken from the title/descriptionof the geneset.")
        ),
        div(
          class = "col-md-4",
          wordcloud_plot_wordtsne_ui(
            ns("wordtsne"),
            height = halfH,
            title = "Word t-SNE",
            info.text = "Keywords that are often found together in title/descriptions are placed close together in the t-SNE. For each keyword we computed enrichment using GSEA on the mean (absolute) enrichment profiles (averaged over all contrasts). Statistically significant gene sets (q<0.05) are colored in red. The sizes of the nodes are proportional to the normalized enrichment score (NES) of the keyword.",
            caption = "T-SNE plot of keywords that were found in the title/description of gene sets. ")
        )
      ),
      div(
        class = "row",
        div(
          class = "col-md-6",
          wordcloud_table_enrichment_ui(
            ns("wordcloud_enrichmentTable"),
            title = "Enrichment table",
            caption = "This table shows the keyword enrichment statistics for the selected contrast. The enrichment is calculated using GSEA for occurance of the keywork in the ordered list of gene set descriptions.",
            info.text = "Keyword enrichment table.",
            height = halfH,            
            width = c("100%", "100%")
          )
        ),
        div(
          class = "col-md-6",
          wordcloud_table_leading_edge_ui(
            ns("wordcloud_leadingEdgeTable"),
            height = halfH,                        
            width = c("100%", "100%"),
            title = "Leading-edge table",
            info.text = "This table contains the input datasets used to create the word cloud.",
            caption = "Keyword leading edge table."
          )
        )
      )
    )
  )
  div(
    boardHeader(title = "Word cloud", info_link = ns("wc_info")),
    tabs
  )
}
