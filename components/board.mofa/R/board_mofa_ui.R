##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

MofaInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    ## data set parameters
    shiny::selectInput(ns("selected_trait"), "Select trait:", choices = NULL),
    shiny::selectInput(ns("selected_module"), "Select module:", choices = NULL),
    shiny::selectInput(ns("selected_factor"), "Select factor:", choices = NULL),
    shiny::selectizeInput(ns("show_types"), "Show datatypes:",
                          choices = NULL, multiple = TRUE),
    shiny::br(),
    bslib::accordion(
      id = ns("data_type_accordion"),
      open = FALSE,
      bslib::accordion_panel(
        "Options",
        icon = icon("cog", lib = "glyphicon"),
        shiny::tagList(
          shiny::selectInput(
            ns("kernel"), "Kernel",
            choices = sort(c("DIABLO", "MOFA", "PCA", "MCIA", "WGCNA", "RGCCA")),
            selected = "MOFA"
          ),
          shiny::conditionalPanel(
            "input.kernel != 'WGCNA'",
            ns = ns,
            shiny::selectInput(ns("numfactors"), "Number of factors",
              choices = c(3, 5, 10, 15, 25),
              selected = 10
            ),
          ),
          shiny::checkboxInput(ns("add_gsets"), "Add gene set layers",
            value = FALSE
          ),
          br(),
          shiny::actionButton(ns("compute"), "Compute!")
        )
      )
    )
  )
}

MofaUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

  #  shiny::div(
  bslib::page_fillable(
    fillable_mobile = FALSE,  # not working here... 
    boardHeader(title = "Multi-Omics Factor Analyis", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),

      ##----------------------------------------------------------------      
      shiny::tabPanel(
        "Overview",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>Multi‐Omics Factor Analysis (MOFA)</b> is a computational, factorization-based framework for multi‐omics data integration. The inferred latent 'factors' (or 'modules') represent the underlying principal axes of heterogeneity across the samples.")),
          bslib::layout_columns(
            col_widths = c(7, 5),
            height = "calc(100vh - 180px)",
            bslib::layout_columns(
              col_widths = c(6, 6),
              mofa_plot_variance_ui(
                ns("factorxview"),
                title = "Variance per factor and type",
                info.text = "Amount of variance explained by each factor in each omic type. A trained MOFA model is used to infer the proportion of variance explained (i.e. the coefficient of determinations (R^2)) by the MOFA factors across the different views. Higher variance suggests stronger effect. In MOFA, 'views' refer to features from non-overlapping set of omic types. MOFA 'factors' are low-dimensional representations of multi-omic data. A factor is a latent variable that captures a source of variation across the integrated data. Each factor captures a different source and dimension of heterogeneity in the integrated data, and thus represents an independent source of variation. Note that the interpretation of factors is analogous to the interpretation of the principal components in PCA. Factors with higher explained variance are typically considered more important for understanding the underlying structure and patterns in a multi-omics dataset. They may correspond to significant biological processes, cellular states, or experimental conditions that have a broader impact across multiple data modalities.",
                info.references = list(list("Argelaguet R, Velten B, Arnol D, Dietrich S, Zenz T, Marioni JC, Buettner F, Huber W, Stegle O (2018). “Multi‐Omics Factor Analysis — a framework for unsupervised integration of multi-omics data sets.” Mol Syst Biol.", "https://www.embopress.org/doi/full/10.15252/msb.20178124")),
                caption = "Amount of variance explained by each factor in each omics type.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              mofa_plot_dendrogram_ui(
                ns("dendrogram"),
                title = "Feature clustering",
                caption = "Factor clustering",                
                info.text = "Gene modules are detected as branches of the resulting cluster tree using the dynamic branch cutting approach. Genes inside a given module are summarized with the module eigengene. The module eigengene of a given module is defined as the first principal component of the standardized expression profiles.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),              
              mofa_plot_variance_ui(
                ns("variance_view"),
                title = "Variance per type",
                info.text = "Total amount of variance explained by each view (omic types). Distinct omic types or data modalities often account for different variance observed in the data. This represents an expectation in multi-omics data analyses. A data type explaining more variance compared to another data type may capture either more biologically heterogenous signals, or be affected by technical noise.",
                info.references = list(list("Argelaguet R, Velten B, Arnol D, Dietrich S, Zenz T, Marioni JC, Buettner F, Huber W, Stegle O (2018). “Multi‐Omics Factor Analysis — a framework for unsupervised integration of multi-omics data sets.” Mol Syst Biol.", "https://www.embopress.org/doi/full/10.15252/msb.20178124")),
                caption = "Total amount of variance explained by each view (omic types). ",
                label = "c",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              mofa_plot_variance_ui(
                ns("variance_factor"),
                title = "Variance per factor",
                info.text = "Amount of variance explained by each factor across all views (omic types) combined. Distinct factors capture different dimensions of heterogeneity in the data. When a factor explains more variance compared to another factor, it is interpreted as a more dominant source of heterogeneity across samples, resulting from stronger biological or technical influence on the data.",
                info.references = list(list("Argelaguet R, Velten B, Arnol D, Dietrich S, Zenz T, Marioni JC, Buettner F, Huber W, Stegle O (2018). “Multi‐Omics Factor Analysis — a framework for unsupervised integration of multi-omics data sets.” Mol Syst Biol.", "https://www.embopress.org/doi/full/10.15252/msb.20178124")),
                caption = "Amount of variance explained by each factor across all views (omic types) combined.",
                label = "c",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              )
            ),
            mofa_plot_moduleheatmap_ui(
              ns("integrated_heatmap"),
              title = "Multi-omics heatmap",
              info.text = "...",
              caption = "Integrated Multi-omics heatmap.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )          
        )
      ),


      ##----------------------------------------------------------------      
      shiny::tabPanel(
        "Response",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>Factor response analysis.</b> Associations of factors with our trait of interest are quantified by the correlation between factor and trait vectors. Factors with high (absolute) factor-trait correlation show large differences between phenotype conditions.")),
          bslib::layout_columns(
            height = "calc(100vh - 180px)",
            col_widths = bslib::breakpoints(
              xxxl = c(12, 12),
              lg = c(6,6),
              sm = c(12, 12)               
            ),
            bslib::layout_columns(
              col_widths = bslib::breakpoints(
                xxxl = c(5, 3, 4),
                lg = c(12, 6, 6),
                sm = c(12, 12, 12, 12, 12)
              ),
              mofa_plot_factortrait_ui(
                ns("factortrait"),
                title = "Factor-Trait correlation",
                info.text = "Correlation between MOFA factors and traits. Correlation between each MOFA factor (y-axis) and each available variable in your metadata (x-axis). The metadata variables shown correspond to those provided in your uploaded samples.csv file. Covariates may include the phenotype(s) of interest. Correlations are pairwise Pearson's correlation coefficients. The colors in the heatmap indicate the strength and direction of the correlation. The correlation values range from -1 to +1. Stronger, positive correlations will approach darker red. Stronger, negative correlations will approach darker blue. Each MOFA factor captures a source of variation across the integrated data types. Therefore, the heatmap helps identification of dominant sources of variation in the data -which may be driven by biological processes- associated with a trait or condition.",
                caption = "Correlation between factors and traits. Correlation between each MOFA factor (y-axis) and each available variable in your metadata (x-axis). Correlations are pairwise Pearson's correlation coefficients. The colors in the heatmap indicate the strength and direction of the correlation. The correlation values range from -1 to +1.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              mofa_plot_factorcorheatmap_ui(
                ns("factorcorheatmap"),
                title = "Correlation between factors.",
                info.text = "Correlation between MOFA factors. Correlations are pairwise Pearson's correlation coefficients. The colors in the heatmap indicate the strength and direction of the correlation. The correlation values range from -1 to +1. Stronger, positive correlations will approach darker red. Stronger, negative correlations will approach darker blue. Correlated factors might capture similar sources of variation in the data -which might be driven by shared biological processes.",
                caption = "",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              mofa_plot_factorgraph_ui(
                ns("factorgraph"),
                title = "Between-factors graph",
                info.text = "For each module, we also define a quantitative measure of module membership (MM) as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module.",
                caption = "",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              )
            ),
            mofa_plot_boxplots_ui(
              ns("boxplots"),
              title = "Factor scores by response",
              info.text = "Distribution of MOFA factors' score per metadata variable. Scores of each selected MOFA factor stratified by a trait/condition/phenotype of interest. The metadata variables shown correspond to those provided in your uploaded samples.csv file. Differences in score distribution between two conditions might suggest distinct sources of variation in each group captured by the selected MOFA factor. In such cases, each group would provide a different contribution to variation in the data captured by the selected MOFA factor.",
              caption = "Distribution of MOFA factors' score per metadata variable.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),
      
      ##----------------------------------------------------------------
      shiny::tabPanel(
        "Weights",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>MOFA weights.</b> For all features in the module, we compute the gene significance (GS) as the correlation of its gene expression and the phenotype (or trait). We also define the 'module membership' (MM) as the correlation of the gene expression and the module eigengene. We can identify 'driver genes' that have high GS as well as high MM.")),
          bslib::layout_columns(
            height = "calc(100vh - 180px)",
            col_widths = bslib::breakpoints(
              lg = c(6,6,6,3,3),
              sm = c(12, 12, 12, 12)
            ),
            mofa_plot_weights_ui(
              ns("weights"),
              title = "Factor weights (aka loadings)",
              info.text = "",
              caption = "Module-trait analysis identifies modules that are significantly associated with the measured traits by quantifying the association as the correlation of the eigengenes with external traits.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            mofa_plot_moduleheatmap_ui(
              ns("module_heatmap"),
              title = "Factor heatmap",
              info.text = "...",
              caption = "Factor heatmap.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            mofa_table_genetable_ui(
              ns("mofa_genetable"),
              title = "Factor features",
              info.text = "",
              caption = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("100%", "100%")
            ),
            mofa_plot_modulegraph_ui(
              ns("modulegraph"),
              title = "Factor feature graph",
              info.text = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            mofa_plot_centrality_ui(
              ns("centrality"),
              title = "Centrality vs. weight",
              info.text = "...",
              caption = "Module heatmap.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )                          
          )
        )
      ), ## tabPanel
 
      ##----------------------------------------------------------------
      shiny::tabPanel(
        "Enrichment",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>Factor enrichment.</b> To understand the biological function for each factor, we can perform a geneset enrichment analysis using its factor loadings. For a selected geneset, the boxes show the joint enrichment, mixed heatmap, pathwaydiagram and feature table.")),
          bslib::layout_columns(
            height = "calc(100vh - 180px)",
            col_widths = bslib::breakpoints(
              lg = c(4, 4, 4, 8, 4),
              sm = c(12, 12, 12, 12)
            ),
            mofa_plot_enrichment_ui(
              ns("plotenrichment"),
              title = "Geneset enrichment analysis of MOFA factors.",
              info.text = "Geneset enrichment analysis of MOFA factors. Please select an enriched gene set for the selected MOFA factor in the Factor enrichment table below. Note that MOFA factors can be selected from the dropdown menu on the right side of the board. In each enrichment plot, black vertical bars indicate the rank of genes in the selected geneset. The green curve corresponds to the profile of the normalized enrichment score computed from the Kolmogovov-Smirnov test. Upper left shift of the green curve indicates enrichment in the first group. Lower right shift of the green curve indicates enrichment in the second group.",
              caption = "Functional analysis of factor",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),            
            mofa_plot_pathwayheatmap_ui(
              ns("pathwayheatmap"),
              title = "Pathway heatmap",
              info.text = "...",
              caption = "Integrated Multi-omics heatmap.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            mofa_plot_pathbank_ui(
              ns("pathway"),
              title = "Pathway diagram",
              caption = "Pathways that integrate proteomics and metabolomics data types in a single pathway diagram.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("100%", "100%")
            ),                        
            mofa_table_factorenrichment_ui(
              ns("mofa_factorenrichment"),
              title = "Factor enrichment table",
              info.text = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            mofa_table_enrichmentgenes_ui(
              ns("mofa_enrichmentgenes"),
              title = "Enrichment gene table",
              label = "a",
              info.text = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )                          
          )
        )
      ) ## tabPanel
      
    )
  )
}
