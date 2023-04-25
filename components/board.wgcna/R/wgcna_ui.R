##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

WgcnaInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),

    ## data set parameters
    shiny::selectInput(ns("selected_module"), "select module", choices = NULL),
    shiny::actionButton(ns("compute"), "Compute!",
      icon = icon("running"),
      class = "btn-outline-primary"
    ),
    shiny::br(),
    shiny::br(),
    shiny::actionLink(ns("options"), "Options", icon = icon("cog", lib = "glyphicon")),
    shiny::br(), br(), br(),
    shiny::conditionalPanel(
      "input.options % 2 == 1",
      ns = ns,
      shiny::tagList(
        shiny::selectInput(ns("ngenes"), "Number genes:",
          choices = c(500, 1000, 2000, 4000, 8000),
          selected = 1000
        ),
        shiny::selectInput(ns("minmodsize"), "Min. module size",
          choices = c(10, 30, 100, 250),
          selected = 30
        ),
        shiny::selectInput(ns("power"), "Power", c(2, 4, 6, 10), selected = 6),
        shiny::selectInput(ns("deepsplit"), "deepsplit", 0:4, selected = 2),
        shiny::selectInput(ns("cutheight"), "Merge cut height",
          choices = c(0.05, 0.10, 0.25, 0.5, 0.9, 0.999),
          selected = 0.25
        )
      )
    )
  )
}

WgcnaUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

  shiny::div(
    # height = 750,
    boardHeader(title = "WGCNA", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),
      shiny::tabPanel(
        "WGCNA",
        bs_alert(HTML("<b>WGCNA module detection.</b> <b>(a)</b> Modules are detected as branches of the resulting cluster tree using the dynamic branch cutting approach. Genes inside a given module are summarized with the module eigengene. The module eigengene of a given module is defined as the first principal component of the standardized expression profiles. <b>(b)</b> Scale independence and mean connectivity plots to determine the soft threshold. <b>(c)</b> Topological overlap matrix visualized as heatmap. <b>(d)</b> Dimensionality reduction maps colored by WGCNA module. <b>(e)</b> Graph network of WGCNA modules.")),
        bslib::layout_column_wrap(
          width = 1,
          height = "calc(100vh - 180px)",
          heights_equal = "row",
          bslib::layout_column_wrap(
            width = 1/2,
            height = "35%",
            wgcna_plot_gdendogram_ui(
              ns("geneDendro"),
              caption = "WGCNA gene dendrogram and gene modules",
              info.text = "Gene modules are detected as branches of the resulting cluster tree using the dynamic branch cutting approach. Genes inside a given module are summarized with the module eigengene. The module eigengene of a given module is defined as the first principal component of the standardized expression profiles.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_s_independence_ui(
              ns("topologyPlots"),
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          ),
          bslib::layout_column_wrap(
            width = 1/3,
            height = "65%",            
            wgcna_plot_TOMheatmap_ui(
              ns("TOMplot"),
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_gclustering_ui(
              ns("umap"),
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_module_graph_ui(
              ns("moduleGraph"),
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),
      shiny::tabPanel(
        "Modules",
        bs_alert(HTML("<b>WGCNA functional analysis.</b> <b>(a)</b> Module-trait analysis identifies modules that are significantly associated with the measured clinical traits by quantifying the association as the correlation of the eigengenes with external traits. <b>(b)</b> Partial correlation network of genes most correlated to the eigengene. <b>(c)</b> Module enrichment plot of top most enriched genesets. <b>(d)</b> Table of genes in the selected module. <b>(e)</b> Functional enrichment of the module calculated using Fisher's exact test.")),
        bslib::layout_column_wrap(
          width = 1,
          height = "calc(100vh - 180px)",
          heights_equal = "row",
          bslib::layout_column_wrap(
            width = 1/3,
            height = "60%",
            wgcna_plot_MTrelationships_ui(
              ns("moduleTrait"),
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_correlation_network_ui(
              ns("corGraph"),
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_correlation_network_ui(
              ns("enrichPlot"),
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          ),
          bslib::layout_column_wrap(
            width = 1,
            height = "40%",            
            style = htmltools::css(grid_template_columns = "1fr 2fr"),
            wgcna_table_genes_ui(
              ns("geneTable"),
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_table_enrichment_ui(
              ns("enrichTable"),
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),
      shiny::tabPanel(
        "Eigengenes",
        bs_alert(HTML("<b>WGCNA eigengene analysis.</b> <b>(a)</b> It is often interesting to visualizing the network of eigengenes and study the relationships among the found modules. One can use the eigengenes as represen- tative profiles and quantify module similarity by eigengene correlation. <b>(b)</b> For each module, we also define a quantitative measure of 'module membership' (MM) as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes to every module.")),
        bslib::layout_column_wrap(
          width = 1,
          height = "calc(100vh - 180px)",
          style = htmltools::css(grid_template_columns = "3fr 9fr"),
          wgcna_plot_eigengene_clustering_ui(
            ns("eigenClustering"),
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          ),
          wgcna_plot_module_membership_ui(
            ns("eigenCorrelation"),
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      ),
      shiny::tabPanel(
        "Intramodular",
        bs_alert(HTML("<b>WGCNA intramodular analysis.</b> We quantify associations of individual genes with our trait of interest (weight) by defining Gene Significance GS as (the absolute value of) the correlation between the gene and the trait. For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules.")
        ),
        bslib::layout_column_wrap(
          width = 1,
          height = "calc(100vh - 180px)",
          style = htmltools::css(grid_template_columns = "3fr 9fr"),
          wgcna_plot_heatmap_membership_ui(
            ns("intraHeatmap"),
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          ),
          wgcna_plot_membership_v_trait_ui(
            ns("intraScatter"),
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )      
        )
      )
    )
  )
}
