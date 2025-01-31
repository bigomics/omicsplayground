##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' DataView module UI Input function
#'
#' @description A shiny Module. Renders the input parts (sidebar contents) for the module.
#'
#' @param id Internal parameters for {shiny}.
#' #'
#' @export
DataViewInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace

  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    withTooltip(shiny::selectInput(ns("search_gene"), tspan("Gene:"), choices = NULL),
      "Type a gene of interest.",
      placement = "top"
    ),
    withTooltip(
      shiny::selectInput(ns("data_samplefilter"), "Filter samples:",
        choices = NULL, multiple = TRUE
      ),
      "Filter the relevant samples for the analysis.",
      placement = "top"
    ),
    withTooltip(shiny::selectInput(ns("data_groupby"), "Group by:", choices = NULL),
      "Select phenotype for grouping the samples.",
      placement = "top"
    ),
    shiny::br(), shiny::br(),
      bslib::accordion(
        id = ns("data_type_accordion"),
        open = FALSE,
        bslib::accordion_panel(
          "Options",
          icon = icon("cog", lib = "glyphicon"),
          withTooltip(
            shiny::radioButtons(
              ns("data_type"), "Scale:",
              choiceNames = c("linear", "log2"), 
              choiceValues = c("counts", "log2"),
              selected = "log2",
              inline = TRUE
            ),
            "Choose an input data type for the analysis.",
            placement = "bottom"
          )
        )
      )
  )
}

#' DataView module UI output function
#'
#' @description Renders the output part for the module as tabsetPanel object
#'
#' @param id Internal parameters for {shiny}.
#' #'
#' @export
DataViewUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  imgH <- c(330, 600) ## heights for small and fullscreen image
  imgH <- c("40vh", "70vh") ## heights for small and fullscreen image
  fullH <- "calc(100vh - 180px)"
  fullH.css <- "height: calc(100vh - 180px);"

  tabs <- shiny::tabsetPanel(
    id = ns("tabs"),
    # Overview tab #####
    shiny::tabPanel(
      "Overview",
      bslib::layout_columns(
        col_widths = 12,
        height = fullH,
        bs_alert("The Overview panel displays data for a selected gene. The 'gene info' box provides more information about the gene and hyperlinks to external databases. The upper plots show the expression level, average expression ranking, and distribution of expression among the samples. The remaining plots, display the most correlated genes and expression in the GTEX tissue database."),
        bslib::layout_columns(
          height = "100%",
          col_widths = c(2, 10),
          dataview_module_geneinfo_ui(
            id = ns("geneinfo"),
            title = "Gene info",
            info.text = "Information about the selected gene and its function from public databases. For more information, follow the hyperlinks to public databases.",
            info.methods = "The available public databases are OMIM [1], KEGG [2] and GO [3]. The R/bioconductor package used to query them is org.Hs.eg.db [4].",
            info.references = list(
              list(
                "Amberger, J. S., Bocchini, C. A., Schiettecatte, F., Scott, A. F., & Hamosh, A. (2015). OMIM. org: Online Mendelian Inheritance in Man (OMIM®), an online catalog of human genes and genetic disorders. Nucleic acids research, 43(D1), D789-D798.",
                "https://doi.org/10.1093/nar/gku1205"
              ),
              list(
                "Kanehisa, M., & Goto, S. (2000). KEGG: kyoto encyclopedia of genes and genomes. Nucleic acids research, 28(1), 27-30.",
                "https://doi.org/10.1093/nar/28.1.27"
              ),
              list(
                "Ashburner, M., Ball, C. A., Blake, J. A., Botstein, D., Butler, H., Cherry, J. M., ... & Sherlock, G. (2000). Gene ontology: tool for the unification of biology. Nature genetics, 25(1), 25-29.",
                "https://doi.org/10.1038/75556"
              ),
              list(
                "Carlson M (2019). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.8.2.",
                "https://bioconductor.org/packages/org.Hs.eg.db/"
              )
            ),
            caption = "Information about the selected gene and its function from public databases.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          ),
          bslib::layout_columns(
            col_widths = 12,
            bslib::layout_columns(
              col_widths = c(4, 4, 4),
              dataview_plot_expression_ui(
                id = ns("expressionplot"),
                title = "Gene expression",
                info.text = "Expression of the selected gene by sample. Samples can be grouped by phenotype using the {Group by} setting. Also, the type of expression can be selected under {Data type} under Settings > Options.",
                caption = "Barplot of abundance or expression of grouped samples (or cells) for the gene selected in the Search gene.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                label = "a"
              ),
              dataview_plot_averagerank_ui(
                ns("averagerankplot"),
                label = "b",
                title = "Average rank",
                info.text = "Rank of the selected gene by decreasing expression. The type of expression can be selected under {Data type} under Settings > Options.",
                caption = "Ranking of the selected gene by decreasing average expression.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              dataview_plot_tsne_ui(
                ns("tsneplot"),
                label = "c",
                title = "t-SNE clustering",
                info.text = "T-SNE clustering of samples (or cells) colored by the expression of the gene selected.",
                info.methods = "Clustering of the samples using the t-distributed stochastic neighbor embedding (t-SNE) method. Performed using the Rtsne R package [1]. The dark blue color represents an over-expression of the selected gene across samples (or cells).",
                info.references = list(
                  list(
                    "Krijthe JH (2015). Rtsne: T-Distributed Stochastic Neighbor Embedding using Barnes-Hut Implementation. R package version 0.17",
                    "https://doi.org/10.32614/CRAN.package.Rtsne"
                  )
                ),
                info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
                caption = "t-SNE of samples colored by expression of selected gene.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              )
            ),
            bslib::layout_columns(
              col_widths = c(7, 5),
              dataview_plot_correlation_ui(
                ns("correlationplot"),
                label = "d",
                title = "Top correlated genes",
                info.text = "Barplot of the top positively and negatively correlated genes with the selected gene.",
                info.methods = "Correlation between selected gene and top positively and negatively correlated genes performed using Pearson correlation method from the R core stats package. Colors are from absolute expression levels of genes, where the low and high expressions range between the light and dark colors, respectively.",
                info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#correlation-analyses",
                caption = "Barplot of the top positively and negatively correlated genes with the selected gene. Darker color corresponds to higher expression of the gene.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              dataview_plot_tissue_ui(
                ns("tissueplot"),
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%"),
                label = "e",
                title = "Tissue expression (GTEX)",
                info.text = "Top 15 expressing tissues for the selected gene in the tissue expression GTEx database. Colors represent tissue clusters. If species other than human, the human ortholog is used to query the GTEX database.",
                info.methods = "Tissue clusters computed using unsupervised clustering. Perfomed using the hierarchical clustering method from the R core stats package, then queried to the GTEx database [1].",
                info.references = list(
                  list(
                    "GTEx Consortium. (2020). The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science, 369(6509), 1318-1330.",
                    "https://doi.org/10.1126/science.aaz1776"
                  )
                ),
                caption = paste("Top 15 expressing tissues for the selected gene in the tissue expression GTEx database. Colors represent tissue clusters. If species other than human, we use the human ortholog to query the GTEX database.")
              )
            )
          )
        )
      )
    ),

    # QC tab #####
    shiny::tabPanel(
      "Sample QC",
      bslib::layout_columns(
        col_widths = 12,
        height = fullH,
        bs_alert("The Sample QC tab provides an overview of several sample-centric quality control metrics. In this QC tab, the total counts per sample and their distribution among the samples are displayed. This is most useful to check the technical quality of the dataset, such as proportion of ribosomal genes."),
        bslib::layout_columns(
          col_widths = c(4, 4, 4),
          dataview_plot_totalcounts_ui(
            ns("counts_total"),
            label = "a",
            title = "Total counts / Number of detected features",
            info.text = "Average total counts by sample. Samples can be grouped using the {Group by} setting.",
            caption = "Barplot of the total counts or number of detected features for each sample (or group).",
            ## caption = "Barplot of the average counts for each group.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          ),
          dataview_plot_boxplot_ui(
            ns("counts_boxplot"),
            title = "Counts boxplots",
            info.text = "Boxplot of counts by sample. Samples can be grouped using the {Group by} setting.",
            caption = "Distribution of total counts per sample/group. The center horizontal bar correspond to the median.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            label = "b"
          ),
          dataview_plot_histogram_ui(
            ns("counts_histplot"),
            title = "Density distribution of counts",
            info.text = "Histogram of the density distribution of total counts per sample. Samples can be grouped using the {Group by} setting.",
            caption = "Density distribution of total counts per sample/group",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%"),
            label = "c"
          )
        ),
        bslib::layout_columns(
          col_widths = c(5, 7),
          dataview_plot_genetypes_ui(
            ns("counts_genetypes"),
            title = "Proportion of major gene types",
            info.text = "Abundance of genetypes on the loaded data. Genetypes can be ribosomal protein genes, kinases or RNA binding motifs, etc. Samples can be grouped using the {Group by} setting.",
            caption = "Barplot showing the proportion of major gene types.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%"),
            label = "d"
          ),
          dataview_plot_abundance_ui(
            ns("counts_abundance"),
            title = "Proportion of major gene types per sample/group",
            info.text = "Barplot showing the percentage of counts in terms of major gene types such as ribosomal protein genes, kinases or RNA binding motifs for each group. Samples can be grouped using the {Group by} setting.",
            caption = "Barplot showing the proportion of counts of major gene types in samples or groups.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            label = "e",
            width = c("auto", "100%")
          )
        )
      )
    ),
    # counts table tab #####
    shiny::tabPanel(
      "Data table",
      bslib::layout_columns(
        col_widths = 12,
        row_heights = list("auto", 1),
        height = fullH,
        bs_alert("In Data table panel, the exact expression values across the samples can be looked up, where genes are ordered by the correlation with respect to the selected gene. Gene-wise average expression of a phenotype sample grouping is also presented in this table."),
        dataview_table_rawdata_ui(
          ns("rawdatatable"),
          title = "Gene expression table",
          info.text = "The column 'rho' reports the correlation with the gene selected in 'Search gene' in the left side bar. If the data type selected is counts, the geometric mean is calculated. The SD column reports the standard deviation of expression across samples (or cells).",
          caption = "The table shows the gene expression values per sample, or average expression values across the groups.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("100%", "100%")
        )
      )
    ),
    # Sample information #####
    shiny::tabPanel(
      "Sample information",
      bslib::layout_columns(
        height = fullH,
        col_widths = 12,
        row_heights = list("auto", 1.33, 1),
        bs_alert("In the Sample information panel, more complete information about samples can be found."),
        bslib::layout_columns(
          width = 6,
          dataview_plot_phenoheatmap_ui(
            ns("phenoheatmap"),
            title = "Phenotype clustering",
            info.text = "Clustering of the available phenotypes across all samples. Samples can be filtered using the {Filter samples} setting.",
            info.methods = "Column ordering performed using unsupervised clustering on a one-hot encoded matrix. Perfomed using the hierarchical clustering method from the R core stats package.",
            caption = "Clustered heatmap of sample information (i.e. phenotype data)",
            height = c("50%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%"),
            label = "a"
          ),
          dataview_plot_phenoassociation_ui(
            ns("phenoassociation"),
            height = c("50%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%"),
            label = "b",
            title = "Phenotype correlation",
            info.text = "Correlation matrix of the phenotypes.",
            info.methods = "Correlation for categorical vs categorical phenotypes performed using fisher test, categorical vs continuous phenotypes performed using Krusall-Wallace and continuous vs continuous phenotypes performed using pearson correlation (all three done using the core R stats package). The size of the dots correspond to the absolute correlation between two phenotype conditions. Red corresponds to positive correlation, blue corresponds to negative correlation.",
            info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#correlation-analyses",
            caption = "Clustered heatmap of phenotype correlation."
          )
        ),
        bslib::layout_columns(
          col_widths = 12,
          dataview_table_samples_ui(
            ns("sampletable"),
            height = c("50%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%"),
            title = "Sample information",
            info.text = "Phenotype variables starting with a 'dot' (e.g. '.cell cycle' and '.gender' ) have been estimated from the data.",
            caption = "Phenotype information about the samples."
          )
        )
      )
    ),

    # contrasts tab #####
    shiny::tabPanel(
      "Contrasts",
      bslib::layout_columns(
        height = fullH,
        col_widths = 12,
        row_heights = list("auto", 1),
        bs_alert("The Contrasts panel shows information about the phenotype comparisons."),
        dataview_table_contrasts_ui(
          ns("contrastTable"),
          title = "Contrast table",
          info.text = "Here, you can check which samples belong to which groups for the different comparisons. Non-zero entries '+1' and '-1' correspond to the group of interest and control group, respectively. Zero or empty entries denote samples not use for that comparison.",
          caption = "Table summarizing the contrasts of all comparisons.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    )
  )

  div(
    boardHeader(title = "Data View", info_link = ns("board_info")),
    tabs
  )
}
