##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

DrugConnectivityInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    withTooltip(shiny::selectInput(ns("dsea_contrast"), "Contrast:", choices = NULL),
      "Select the contrast corresponding to the comparison of interest.",
      placement = "top"
    ),
    withTooltip(shiny::selectInput(ns("dsea_method"), "Analysis type:", choices = ""),
      "Select type of drug enrichment analysis: activity or sensitivity (if available).",
      placement = "top"
    ),
    shiny::hr(),
    withTooltip(
      shiny::checkboxInput(
        ns("dseatable_filter"),
        "only annotated drugs",
        FALSE
      ),
      "Show only annotated drugs."
    )
  )
}

DrugConnectivityUI <- function(id) {
  ns <- shiny::NS(id)

  fullH <- "calc(100vh - 180px)"
  halfH <- "calc(50vh  - 98px)"

  panel1 <- shiny::tabPanel(
    "Drug enrichment",
    bslib::layout_columns(
      col_widths = c(9, 3),
      height = fullH,
      bslib::layout_columns(
        col_widths = 12,
        bslib::layout_columns(
          col_widths = c(6, 6),
          drugconnectivity_plot_enplots_ui(
            id = ns("dsea_enplots"),
            title = "Drug connectivity",
            info.text = "Not available for this plot",
            caption = "Drug connectivity correlates your signature with known drug profiles from the L1000 database, and shows similar and opposite profiles by running the GSEA algorithm on the drug profile correlation space.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%"),
            label = "a"
          ),
          drugconnectivity_plot_moa_ui(
            id = ns("dsea_moaplot"),
            title = "Mechanism of action",
            info.text = "This plot visualizes the mechanism of action (MOA) across the enriched drug profiles. On the vertical axis, the GSEA normalized enrichment score of the MOA class or gene target is plotted. You can switch to visualize between MOA class or target gene.",
            caption = "Mechanism of action (MOA) plot indicating the most correlated drug MOAs for a selected contrast.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            label = "c"
          )
        ),
        drugconnectivity_table_dsea_ui(
          ns("dsea_table"),
          title = "Enrichment table",
          info.text = "The platform correlates your signature with known drug or single gene alteration profiles from the selected database, and shows similar and opposite profiles by running the GSEA algorithm on the drug or gene alteration profile correlation space. Interpretation of the correlation is similar to standard GSEA plots.",
          caption = "GSEA-like plots showing the correlation of various drug or single gene alteration expression profiles with the selected contrast signature.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("100%", "100%")
        )
      ),
      drugconnectivity_plot_actmap_ui(
        ns("dsea_actmap"),
        title = "Activation matrix",
        info.text = "The Activation Matrix visualizes the activation of drug activation enrichment across the conditions. The size of the circles correspond to their relative activation, and are colored according to their upregulation (red) or downregulation (blue) in the contrast profile.",
        caption = "Activation Matrix visualising drug or single gene alteration profile correlations in all the available contrasts.",
        height = c("100%", TABLE_HEIGHT_MODAL),
        width = c("100%", "100%"),
        label = "d"
      )
    )
  ) ## end of tabPanel

  panel2 <- shiny::tabPanel(
    "Connectivity map (beta)",
    bslib::layout_columns(
      col_widths = c(5, 7),
      bslib::layout_columns(
        col_widths = 12,
        drugconnectivity_plot_cmap_enplot_ui(
          id = ns("cmap_enplot"),
          title = "Enrichment plot",
          info.text = "Not available for this plot.",
          caption = "Enrichment of the selected drug perturbation profile with your selected signature.",
          label = "a",
          height = c(halfH, TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        ),
        drugconnectivity_table_cmap_ui(
          id = ns("cmap_table"),
          title = "Connectivity table",
          info.text = "Enrichment is calculated by correlating your signature with known drug or single gene alteration profiles from teh selected database. Because the databases have multiple perturbation experiments for a single drug or gene, drugs or genes are scored by running the GSEA algorithm on the contrast-drug/gene profile correlation space. In this way, we obtain a single score for multiple profiles of a single drug or gene alteration.",
          caption = "Enrichment table showing the normalised enrichment score and p-values of a selected contrast signature against drug or gene alteration profiles.",
          height = c(halfH, TABLE_HEIGHT_MODAL),
          width = c("100%", "100%")
        )
      ),
      bslib::layout_columns(
        col_widths = 12,
        drugconnectivity_plot_cmap_dsea_ui(
          id = ns("cmap_dsea"),
          title = "Connectivity Map",
          info.text = " The platform correlates your signature with known drug profiles from the L1000 database, and shows similar and opposite profiles by running the GSEA algorithm on the drug profile correlation space. The color corresponds to the rank correlation between the drug signatures and your selected contrast.",
          caption = "Plot showing the top signatures as UMAP. Each point is one L1000 experiment.",
          label = "c",
          height = c(fullH, TABLE_HEIGHT_MODAL)
        )
      )
    )
  )


  div(
    boardHeader(title = "Drug Connectivity", info_link = ns("dsea_info")),
    panel1 <- shiny::tabsetPanel(
      id = ns("tabs"),
      panel1,
      panel2
    )
  )
}
