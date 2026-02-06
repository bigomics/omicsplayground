##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

WordCloudInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
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
  halfH <- c("calc(50vh - 70px)", TABLE_HEIGHT_MODAL)

  ns <- shiny::NS(id) ## namespace
  tabs <- shiny::tabsetPanel(
    id = ns("tabs"),
    shiny::tabPanel(
      "Word Cloud",
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        wordcloud_plot_enrichment_ui(
          ns("gseaplots"),
          title = "Enrichment plots",
          info.text = "Enrichment plots for all available contrasts for genesets that contain the keyword selected on the Enrichment table.",
          info.methods = "Keyword enrichment is computed by running GSEA [1] on the enrichment score profile for all contrasts. The test is defined as the collection of genesets that contain the keyword in the title/description. Black vertical bars indicate the position of gene sets that contains the keyword in the ranked list of enrichment scores. The curve in green corresponds to the 'running statistic' of the keyword enrichment score. The more the green ES curve is shifted to the upper left of the graph, the more the keyword is enriched in the first group. Conversely, a shift of the green ES curve to the lower right, corresponds to keyword enrichment in the second group.",
          info.references = list(
            list(
              "Shi, J., & Walker, M. G. (2007). Gene set enrichment analysis (GSEA) for interpreting gene expression profiles. Current Bioinformatics, 2(2), 133-137.",
              "https://doi.org/10.2174/157489307780618231"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#statistical-testing",
          caption = "Keyword enrichment plots for all available contrasts.",
          height = halfH
        ),
        wordcloud_plot_wordcloud_ui(
          ns("wordcloud"),
          height = halfH,
          title = "Word cloud",
          info.text = "Word cloud of the most enriched keywords for the data set. Word can be excluded using the {Exclude words} plot setting; additionally the color palette of the plot can be selected on the {Colors} plot setting.",
          info.methods = "Keyword enrichment is computed by running GSEA [1] on the mean (squared) enrichment profile (averaged over all contrasts). For each keyword, we defined the 'keyword set' as the collection of genesets that contain that keyword in the title/description. The sizes of the words are relative to the normalized enrichment score (NES) from the GSEA computation.",
          info.references = list(
            list(
              "Shi, J., & Walker, M. G. (2007). Gene set enrichment analysis (GSEA) for interpreting gene expression profiles. Current Bioinformatics, 2(2), 133-137.",
              "https://doi.org/10.2174/157489307780618231"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#statistical-testing",
          caption = "Word cloud of the most enriched keywords for the data set. Words are taken from the title/descriptionof the geneset."
        ),
        wordcloud_plot_wordtsne_ui(
          ns("wordtsne"),
          height = halfH,
          title = "Word t-SNE",
          info.text = "Clustering plot of the keywords of the genesets. The clustering method can be selected under the {Clustering algorithm} plot setting.",
          info.methods = "For each keyword enrichment is computed using GSEA [1] on the mean (absolute) enrichment profiles (averaged over all contrasts). Statistically significant gene sets (q<0.05) are colored in blue. The sizes of the nodes are proportional to the normalized enrichment score (NES) of the keyword. Keywords that are often found together in title/descriptions are placed close together on the plot.",
          info.references = list(
            list(
              "Shi, J., & Walker, M. G. (2007). Gene set enrichment analysis (GSEA) for interpreting gene expression profiles. Current Bioinformatics, 2(2), 133-137.",
              "https://doi.org/10.2174/157489307780618231"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#statistical-testing",
          caption = "Clustering plot of keywords that were found in the title/description of gene sets."
        )
      ),
      bslib::layout_columns(
        col_widths = c(6, 6),
        wordcloud_table_enrichment_ui(
          ns("wordcloud_enrichmentTable"),
          title = "Enrichment table",
          caption = "This table shows the keyword enrichment statistics for the selected contrast. The enrichment is calculated using GSEA for occurance of the keywork in the ordered list of gene set descriptions.",
          info.text = "Keyword enrichment table.",
          height = halfH,
          width = c("100%", "100%")
        ),
        wordcloud_table_leading_edge_ui(
          ns("wordcloud_leadingEdgeTable"),
          height = halfH,
          width = c("100%", "100%"),
          title = "Leading-edge table",
          info.text = "This table contains the input datasets used to create the word cloud.",
          caption = "Keyword leading edge table."
        )
      )
    ),
    shiny::tabPanel(
      "AI Summary",
      bslib::layout_columns(
        col_widths = 12,
        height = "calc(100vh - 181px)",
        wordcloud_ai_summary_ui(
          ns("wordcloud_ai_summary"),
          title = "AI Keyword Summary",
          info.text = "AI-generated interpretation of keyword enrichment results for the selected contrast.",
          caption = "AI-generated keyword enrichment summary.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    )
  )
  div(
    boardHeader(title = "Word cloud", info_link = ns("wc_info")),
    tabs
  )
}
