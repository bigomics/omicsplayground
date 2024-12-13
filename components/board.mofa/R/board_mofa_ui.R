##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

MofaInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),

    ## data set parameters
    shiny::selectInput(ns("selected_factor"), "Select factor", choices = NULL),
    shiny::selectizeInput(ns("show_types"), "Show datatypes",
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
          choices = c("MOFA","PCA","DIABLO","MCIA","wmfcna"),
          selected = "PCA"                 
        ),
        shiny::selectInput(ns("numfactors"), "Number of factors",
          choices = c(3,5,10,15,25),
          selected = 10
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
                info.text = "WGCNA Topological Overlap Matrix (TOM) heatmap.",
                caption = "Topological overlap matrix visualized as heatmap.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              mofa_plot_dendrogram_ui(
                ns("dendrogram"),
                title = "Feature clustering",
                caption = "MOFA factor heatmap",                
                info.text = "Gene modules are detected as branches of the resulting cluster tree using the dynamic branch cutting approach. Genes inside a given module are summarized with the module eigengene. The module eigengene of a given module is defined as the first principal component of the standardized expression profiles.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),              
              mofa_plot_variance_ui(
                ns("variance_view"),
                title = "Variance per type",
                info.text = "WGCNA Topological Overlap Matrix (TOM) heatmap.",
                caption = "Topological overlap matrix visualized as heatmap.",
                label = "c",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              mofa_plot_variance_ui(
                ns("variance_factor"),
                title = "Variance per factor",
                info.text = "WGCNA Topological Overlap Matrix (TOM) heatmap.",
                caption = "Topological overlap matrix visualized as heatmap.",
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
          bs_alert(HTML("<b>MOFA factor response.</b> We quantify associations of factors with our trait of interest (weight) by the correlation between the factor and the trait. For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules.")),
          bslib::layout_columns(
            col_widths = breakpoints(
              xxxl = c(12, 12),
              lg = c(7,5),
              sm = c(12, 12)               
            ),
            bslib::layout_columns(
              col_widths = breakpoints(
                xxxl = c(5, 3, 4),
                lg = c(12, 5, 7),
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
              mofa_plot_moduleheatmap_ui(
                ns("module_heatmap"),
                title = "Factor heatmap",
                info.text = "...",
                caption = "Factor heatmap.",
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
        "Factor",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>MOFA weights.</b> <b>(a)</b>")),
          bslib::layout_columns(
            col_widths = breakpoints(
              lg = c(7, 5, 4, 4, 4),
              sm = c(12, 12, 12, 12, 12)
            ),
            mofa_plot_weights_ui(
              ns("weights"),
              title = "Factor loading",
              info.text = "",
              caption = "Module-trait analysis identifies modules that are significantly associated with the measured traits by quantifying the association as the correlation of the eigengenes with external traits.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            mofa_plot_enrichment_ui(
              ns("enrichment"),
              title = "Functional enrichment analysis",
              info.text = "",
              caption = "Functional analysis of factor",
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
            ),
            mofa_plot_modulegraph_ui(
              ns("modulegraph"),
              title = "Features graph",
              info.text = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),              
            mofa_plot_centrality_ui(
              ns("centrality"),
              title = "Feature importance",
              info.text = "...",
              caption = "Module heatmap.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )              
          )
        )
      ) ## tabPanel
      
    )
  )
}
