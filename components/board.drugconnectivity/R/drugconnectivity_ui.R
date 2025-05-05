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

  fullH <- "calc(100vh - 181px)"
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
            info.text = "GSEA enrichment-like plots of drug profiles correlated with feature signature for the selected {Contrast}.",
            info.methods = "The log2 fold-change (log2FC) matrix is extracted for all computed contrasts. Common genes between the log2FC matrix and L1000 database are identified. If less than 20 genes are shared, the analysis is not conducted due to lack of statistical power. Drug meta sets are defined from available drugs by assigning each drug into principal drug classes. Only drug classes containing at least 10 distinct drug reports are retained for analyses. For expression data and L1000 drug matrix the rank of genes is calculated for each sample and drug, respectively, using the colRanks function from the matrixStats R package, with 'average' as method to treat ties. Rank pairwise Pearson's correlation (R) is then computed between the two rank matrices. Small, random gaussian noise is added to the correlation matrix. A sparse, binary (absence/presence) drug model matrix is created and rank correlation between the rank correlation matrix R and the drug model matrix is computed to assess potential association with drug profile and drug-response. A p-value and FDR are also computed for each drug. Drug set enrichment is also performed for each contrast using fgsea on the pre-defined drug meta sets and the rank correlation matrix. GSEA normalized enrichment score (NES) along with p and q-values for each drug are computed. For both rank correlation and GSEA analysis, the drugs are ranked (by correlation or NES) to identify the top 1000 matching drugs for each contrast. Further analyses and visualizations are conducted using the rank correlation values (R) and GSEA NES scores.",
            info.references = list(
              list(
                "Subramanian, A., Tamayo, P.,  Mootha, V. K., .... & Mesirov, J. P. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. PNAS, 102(43), 15545-15550.",
                "https://www.pnas.org/doi/10.1073/pnas.0506580102"
              ),
              list(
                "Korotkevich, G., Sukhov, V., Budin, N., .... & Sergushichev, A. (2021). Fast gene set enrichment analysis BioRxiv.",
                "https://www.biorxiv.org/content/10.1101/060012v3"
              )
            ),
            caption = "Drug connectivity correlates your signature with known drug profiles from the L1000 database. Similar and opposite profiles are shown by running the GSEA algorithm on the drug profile correlation space. As rank metrix (y-axis), the rank pairwise Pearson's correlation matrix (R) is used. See 'Methods' for details on how this is calculated.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%"),
            label = ""
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
          info.text = "Enrichment table summarizing the statistical results of GSEA-based drug enrichment analysis. Enrichment is calculated by correlating your signature with known drug profiles from the L1000 database. Because the L1000 has multiple perturbation experiment for a single drug, drugs are scored by running the GSEA algorithm on the contrast-drug profile correlation space. In this way, we obtain a single score for multiple profiles of a single drug. For details on how rank correlation and GSEA NES values are computed, please refer to the information provided for the 'Drug connectivity' plots.",
          caption = "Drug enrichment table. GSEA-like plots showing the correlation of various drug or single gene alteration expression profiles with the selected contrast signature.",
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
