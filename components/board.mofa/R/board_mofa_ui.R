##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

MofaInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),

    ## data set parameters
    shiny::selectInput(ns("selected_factor"), "Select factor:", choices = NULL),
    shiny::selectInput(ns("selected_module"), "Select module:", choices = NULL),    
    shiny::selectizeInput(ns("show_types"), "Show datatypes:",
                          choices=NULL, multiple=TRUE),
    shiny::br(),
    shiny::br(),
    shiny::br(),    
    shiny::actionLink(ns("options"), "Options", icon = icon("cog", lib = "glyphicon")),
    shiny::br(), 
    shiny::conditionalPanel(
      "input.options % 2 == 1",
      ns = ns,
      shiny::tagList(
        shiny::selectInput(
          ns("kernel"), "Kernel",
          choices = c("MOFA","PCA","DIABLO","MCIA","WGCNA"),
          selected = "MOFA"                 
        ),
        shiny::conditionalPanel(
          "input.kernel != 'WGCNA'",
          ns = ns,
          shiny::selectInput(ns("numfactors"), "Number of factors",
            choices = c(3,5,10,15,25),
            selected = 10
          )
        ),
        shiny::checkboxInput(ns("add_gsets"), "Add gene set layers",
          value = FALSE
        ),
        br(),
        shiny::actionButton(ns("compute"),"Compute!")
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
          bs_alert(HTML("<b>Multi‐Omics Factor Analysis (MOFA)</b> is a computational framework for multi‐omics data integration. The inferred latent 'factors' (or 'modules') represent the underlying principal axes of heterogeneity across the samples.")),
          bslib::layout_columns(
            col_widths = c(7, 5),
            bslib::layout_columns(
              col_widths = c(6, 6),
              mofa_plot_variance_ui(
                ns("factorxview"),
                title = "Variance per factor and type",
                info.text = "Amount of variance explained by each factor in each omic type. A trained MOFA model is used to infer the proportion of variance explained (i.e. the coefficient of determinations (R^2)) by the MOFA factors across the different views. Higher variance suggests stronger effect. In MOFA, 'views' refer to features from non-overlapping set of omic types. MOFA 'factors' are low-dimensional representations of multi-omic data. A factor is a latent variable that captures a source of variation across the integrated data. Each factor captures a different source and dimension of heterogeneity in the integrated data, and thus represents an independent source of variation. Note that the interpretation of factors is analogous to the interpretation of the principal components in PCA. Factors with higher explained variance are typically considered more important for understanding the underlying structure and patterns in a multi-omics dataset. They may correspond to significant biological processes, cellular states, or experimental conditions that have a broader impact across multiple data modalities.",
                caption = "Amount of variance explained by each factor in each omics type.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              mofa_plot_dendrogram_ui(
                ns("dendrogram"),
                title = "Feature clustering",
                caption = "Factor heatmap",                
                info.text = "Gene modules are detected as branches of the resulting cluster tree using the dynamic branch cutting approach. Genes inside a given module are summarized with the module eigengene. The module eigengene of a given module is defined as the first principal component of the standardized expression profiles.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),              
              mofa_plot_variance_ui(
                ns("variance_view"),
                title = "Variance per type",
                info.text = "Total amount of variance explained by each view (omic types). Distinct omic types or data modalities often account for different variance observed in the data. This represents an expectation in multi-omics data analyses. A data type explaining more variance compared to another data type may capture either more biologically heterogenous signals, or be affected by technical noise.",
                caption = "Total amount of variance explained by each view (omic types). ",
                label = "c",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              mofa_plot_variance_ui(
                ns("variance_factor"),
                title = "Variance per factor",
                info.text = "Amount of variance explained by each factor across all views (omic types) combined. Distinct factors capture different dimensions of heterogeneity in the data. When a factor explains more variance compared to another factor, it is interpreted as a more dominant source of heterogeneity across samples, resulting from stronger biological or technical influence on the data.",
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
          bs_alert(HTML("<b>Factor response analysis.</b> We quantify associations of factors with our trait of interest (weight) by the correlation between the factor and the trait. For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules.")),
          bslib::layout_columns(
            col_widths = bslib::breakpoints(
              xxxl = c(12, 12),
              lg = c(7,5),
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
                info.text = "",
                caption = "",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              mofa_plot_factorheatmap_ui(
                ns("factor_heatmap"),
                title = "Between-factor correlation",
                info.text = "For each module, we also define a quantitative measure of module membership (MM) as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module.",
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
              info.text = "Partial correlation graph centered on module eigen-gene with top most correlated features. Green edges correspond to positive (partial) correlation, red edges to negative (partial) correlation. Width of the edges is proportional to the correlation strength of the gene pair. The regularized partial correlation matrix is computed using the 'graphical lasso' (Glasso) with BIC model selection.",
              caption = "Module enrichment plot of top most enriched genesets.",
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
          bs_alert(HTML("<b>MOFA weights.</b> <b>(a)</b>")),
          bslib::layout_columns(
            col_widths = bslib::breakpoints(
              lg = c(6,6,6,3,3),
              sm = c(12, 12, 12, 12)
            ),
            mofa_plot_weights_ui(
              ns("weights"),
              title = "Factor weights",
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
          bs_alert(HTML("<b>MOFA weights.</b> <b>(a)</b>")),
          bslib::layout_columns(
            col_widths = bslib::breakpoints(
              lg = c(6, 6, 12),
              sm = c(12, 12, 12, 12)
            ),
            mofa_plot_enrichment_ui(
              ns("enrichmentplot"),
              title = "Factor enrichment",
              info.text = "",
              caption = "Functional analysis of factor",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),            
            mofa_plot_pathbank_ui(
              ns("pathway"),
              title = "Multi-omics pathway",
              caption = "Pathways that integrate proteomics and metabolomics data types in a single pathway diagram.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("100%", "100%")
            ),                        
            mofa_table_enrichment_ui(
              ns("mofa_enrichmenttable"),
              title = "Enrichment table",
              info.text = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )                          
          )
        )
      ), ## tabPanel

      ##----------------------------------------------------------------
      shiny::tabPanel(
        "gsetMOFA",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>Geneset MOFA.</b> Instead of performing factor analysis on the genes, we can use MOFA directly on a (single-samples) geneset expression matrix. In this way we can find 'genesets modules' (GM) and similarly correlated the modules with traits. The GenesetModule-Factor correlation shows how modules correlate with (gene) factors.")),
          bslib::layout_columns(
            col_widths = bslib::breakpoints(
              lg = c(6, 6, 12),
              sm = c(12, 12, 12, 12)
            ),
            mofa_plot_gsetmofa_traitCor_ui(
              ns("gset_traitcor"),
              title = "GenesetModule-Trait correlation",
              info.text = "...",
              caption = "Factor heatmap.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            mofa_plot_gsetmofa_factorCor_ui(
              ns("gset_factorcor"),
              title = "GenesetModule-Factor correlation",
              info.text = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            mofa_table_gsetmofa_ui(
              ns("gsetmofa_table"),
              title = "Geneset MOFA table",
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
