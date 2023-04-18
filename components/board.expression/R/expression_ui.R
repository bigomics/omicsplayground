##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


ExpressionInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    withTooltip(shiny::selectInput(ns("gx_contrast"), "Contrast:", choices = NULL),
      "Select a contrast of interest for the analysis.",
      placement = "top"
    ),
    withTooltip(shiny::selectInput(ns("gx_features"), "Gene family:", choices = NULL, multiple = FALSE),
      "Choose a specific gene family for the analysis.",
      placement = "top"
    ),
    shiny::fillRow(
      flex = c(1, 1),
      withTooltip(shiny::selectInput(ns("gx_fdr"), "FDR", choices = c(1e-9, 1e-6, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, 1), selected = 0.2),
        "Set the false discovery rate (FDR) threshold.",
        placement = "top"
      ),
      withTooltip(
        shiny::selectInput(ns("gx_lfc"), "logFC",
          choices = c(0, 0.1, 0.2, 0.5, 1, 2, 5), selected = 0
        ),
        "Set the logarithmic fold change (logFC) threshold.",
        placement = "top"
      )
    ),
    shiny::br(), br(), br(), br(),
    withTooltip(shiny::actionLink(ns("gx_options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top"
    ),
    shiny::br(), br(),
    shiny::conditionalPanel(
      "input.gx_options % 2 == 1",
      ns = ns,
      shiny::tagList(
        withTooltip(shiny::checkboxInput(ns("gx_showall"), "show all genes", FALSE),
          "Display all genes in the table. Disable filtering of significant genes.",
          placement = "top", options = list(container = "body")
        ),
        withTooltip(
          shiny::checkboxGroupInput(ns("gx_statmethod"), "Statistical methods:",
            choices = NULL, inline = TRUE
          ),
          "Select a method for the statistical test. To increase the statistical reliability of the Omics Playground,
           we perform the DE analysis using commonly accepted methods in the literature, including t-test (standard, 
           Welch), limma (no trend, trend, voom), edgeR (QLF, LRT), and DESeq2 (Wald, LRT), and merge the results.",
          placement = "right", options = list(container = "body")
        )
      )
    )
  )
}

ExpressionUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- "calc(100vh - 120px)"  ## full height of page (minus header)
  halfH <- "height: calc(50vh - 120px);" ## half height of page
  
  tabs1 <- shiny::tabsetPanel(
    id = ns("tabs1"),
    shiny::tabPanel(
      "Plot",
      bslib::layout_column_wrap(
        width = 1/4,
        style = halfH,
        expression_plot_volcano_ui(ns("plots_volcano"),
          label = "a",
          title = "Volcano plot",
          info.text = "A volcano plot of genes for the selected comparison under the Contrast settings plotting fold-change versus significance on the x and y axes, respectively.",
          caption = "Volcano-plot displaying fold-change versus significance.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        ),
        expression_plot_maplot_ui(
          id = ns("plots_maplot"),
          title = "Bland-Altman (MA) plot",
          info.text = "An application of a Bland-Altman (MA) plot of genes for the selected comparison under the Contrast settings plotting mean intensity versus fold-change on the x and y axes, respectively.",
          caption = "MA-plot displaying signal intensity versus fold-change.",
          label = "b",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        ),
        expression_plot_barplot_ui(
          id = ns("plots_barplot"),
          title = "Differential expression",
          info.text = "The top N = {12} differentially (both positively and negatively) expressed gene barplot for the selected comparison under the Contrast settings.",
          caption = "Sorted barplot of the top diffentially expressed genes with largest (absolute) fold-change for selected contrast.",
          label = "c",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        ),
        expression_plot_topfoldchange_ui(
          id = ns("plots_topfoldchange"),
          title = "Gene in contrasts",
          info.text = "The fold change summary barplot across all contrasts for a gene that is selected from the differential expression analysis table under the Table section.",
          caption = "Sorted barplot of the differential expression of the selected gene across all contrasts.",
          label = "d",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    ),
    shiny::tabPanel(
      "Top genes",
      bslib::layout_column_wrap(
        width = 1,
        style = halfH,        
        expression_plot_topgenes_ui(
          id = ns("topgenes"),
          title = "Expression of top differentially expressed genes",
          info.text = "Under the plot Settings, users can scale the abundance levels (counts) or ungroup the samples in the plot from the log scale and ungroup samples settings, respectively.",
          caption = "Expression barplots of the top most differentially expressed genes for the selected contrast.",
          label = "a",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    ),
    shiny::tabPanel(
      "Volcano (all)",
      bslib::layout_column_wrap(
        width = 1,
        style = halfH,        
        expression_plot_volcanoAll_ui(
          id = ns("volcanoAll"),
          title = "Volcano plots for all contrasts",
          info.text = "Under the Volcano (all) tab, the platform simultaneously displays multiple volcano plots for genes across all contrasts. This provides users an overview of the statistics for all comparisons. By comparing multiple volcano plots, the user can immediately see which comparison is statistically weak or strong.",
          caption = "Simultaneous visualisation of volcano plots of genes for all contrasts.",
          label = "a",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    ),
    shiny::tabPanel(
      "Volcano (methods)",
      bslib::layout_column_wrap(
        width = 1,
        style = halfH,                
        expression_plot_volcanoMethods_ui(
          id = ns("volcanoMethods"),
          title = "Volcano plots for all methods",
          info.text = "These plots provide users an overview of the statistics of all methods at the same time. Methods showing better statistical significance will show volcano plots with 'higher wings'.",
          caption = "Simultaneous visualisation of volcano plots of genes by multiple differential expression methods for the selected contrast. ",
          label = "a",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    ) ## end upper tabPanel
  ) ## end of tabs1
  
  tabs2 <- shiny::tabsetPanel(
    id = ns("tabs2"),
    shiny::tabPanel(
      "Table",
      bslib::layout_column_wrap(
        width = 1,
        style = paste(halfH, htmltools::css(grid_template_columns = "2fr 1fr")),
        expression_table_genetable_ui(
          ns("genetable"),
          title = "Differential expression analysis",
          info.text = "The table shows the results of the statistical tests. To increase the statistical reliability of the Omics Playground, we perform the DE analysis using four commonly accepted methods in the literature, namely, T-test (standard, Welch), limma (no trend, trend, voom), edgeR (QLF, LRT), and DESeq2 (Wald, LRT), and merge the results. For a selected comparison under the Contrast setting, the results of the selected methods are combined and reported under the table, where meta.q for a gene represents the highest q value among the methods and the number of stars for a gene indicate how many methods identified significant q values (q < 0.05). The table is interactive (scrollable, clickable); users can sort genes by logFC, meta.q, or average expression in either conditions. Users can filter top N = {10} differently expressed genes in the table by clicking the top 10 genes from the table Settings.",
          caption = "Table showing the significant results of the differential expression analysis on the selected contrast.",
          width = c("100%", "100%"),
          height = c("100%", TABLE_HEIGHT_MODAL)
        ),
        expression_table_gsettable_ui(
          ns("gsettable"),
          title = "Gene sets with gene",
          info.text = "By clicking on a gene in the Table I, it is possible to see which genesets contain that gene in this table, and check the differential expression status in other comparisons from the Gene in contrasts plot under the Plots tab.",
          caption = "Table indicating all the gene sets that contain the gene highlighted in the differential expression table.",
          width = c("100%", "100%"),
          height = c("100%", TABLE_HEIGHT_MODAL)
        )
      ) ## end of layout
    ),
    shiny::tabPanel(
      "Foldchange (all)",
      bslib::layout_column_wrap(
        width = 1,
        style = halfH,                  
        expression_table_fctable_ui(
          ns("fctable"),
          title = "Gene fold changes for all contrasts",
          info.text = "Differential expression (fold-change) across all contrasts. The column `rms.FC` corresponds to the root-mean-square fold-change across all contrasts.",
          caption = "The Foldchange (all) table reports the gene fold changes for all contrasts in the selected dataset.",
          width = c("100%", "100%"),
          height = c("100%", TABLE_HEIGHT_MODAL)
        )
      )
    ),
    shiny::tabPanel(
      "FDR table",
      bslib::layout_column_wrap(
        width = 1,
        style = halfH,          
        expression_table_FDRtable_ui(
          ns("FDRtable"),
          title = "Number of significant genes",
          info.text = "Use this table to quickly see which methods are more sensitive. The left part of the table (in blue) correspond to the number of significant down-regulated genes, the right part (in red) correspond to the number of significant overexpressed genes.",
          caption = "This table reports the number of significant genes at different FDR thresholds for all contrasts and methods.",
          width = c("100%", "100%"),
          height = c("100%", TABLE_HEIGHT_MODAL)
        )
      )
    )      
  ) ## end tabs2 = bottom tabsetPanel

  div(
    boardHeader(title = "Differential expression", info_link = ns("gx_info")),
    bslib::layout_column_wrap(
      width = 1,
      height = fullH,
      gap = '0px',
      tabs1,
      tabs2
    )
  )


}
