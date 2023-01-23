##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


ExpressionInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    withTooltip(
      shiny::actionLink(ns("gx_info"), "Tutorial", icon = shiny::icon("youtube")),
      "Show more information about this module."
    ),
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
                        we perform the DE analysis using commonly accepted methods in the literature, including t-test (standard, Welch),
                        limma (no trend, trend, voom), edgeR (QLF, LRT), and DESeq2 (Wald, LRT), and merge the results.",
          placement = "right", options = list(container = "body")
        )
      )
    )
  )
}

ExpressionUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 800 ## full height of page
  rowH <- 340 ## full height of page
  imgH <- 340 ## height of images

  tagList(
    div(
      style = "max-height:50vh;",
      shiny::tabsetPanel(
        id = ns("tabs1"),
        shiny::tabPanel(
          "Plot",
          div(
            class = "row",
            div(
              class = "col-md-3",
              expression_plot_volcano_ui(ns("plots_volcano"),
                label = "A",
                height = c(imgH, imgH),
                width = c("auto", imgH)
              ),
            ),
            div(
              class = "col-md-3",
              expression_plot_maplot_ui(
                id = ns("plots_maplot"),
                label = "B",
                height = c(imgH, imgH),
                width = c("auto", imgH)
              ),
            ),
            div(
              class = "col-md-3",
              expression_plot_boxplot_ui(
                id = ns("plots_boxplot"),
                label = "C",
                height = c(imgH, imgH),
                width = c("auto", imgH)
              ),
            ),
            div(
              class = "col-md-3",
              expression_plot_topfoldchange_ui(
                id = "plots_topfoldchange",
                label = "D",
                height = c(imgH, imgH),
                width = c("auto", imgH)
              ),
            )
          ),
          tags$div(
            HTML("<b>Expression plots</b> associated with the selected contrast. <b>(a)</b> Volcano-plot plotting fold-change versuson
                            significance the x and y axes, respectively. <b>(b)</b> MA-plot plotting signal intensity versus fold-change on the x and y axes,
                            respectively. <b>(c)</b> Sorted barplot of the top diffentially expressed genes with largest (absolute) fold-change
                            for selected contrast. <b>(d)</b> Sorted barplot of the differential expression of the selected gene across all contrasts.")
          )
        ),
        shiny::tabPanel(
          "Top genes",
          expression_plot_topgenes_ui(
            id = "topgenes",
            label = "A",
            height = c(0.45 * fullH, 700), # c(imgH,420)
            width = c("auto", 1200)
          ), # c('auto',1600)

          shiny::br(),
          tags$div(
            HTML("<b>Top differentially expressed genes.</b> Expression barplots of the top most differentially
                            (both positively and negatively) expressed genes for the selected contrast.")
          )
        ),
        shiny::tabPanel(
          "Volcano (all)",
          plotWidget(ns("volcanoAll")),
          shiny::br(),
          tags$div(
            HTML("<b>Volcano plot for all contrasts.</b> Simultaneous visualisation of volcano
                        plots of genes for all contrasts. Experimental contrasts with better statistical significance will
                        show volcano plots with 'higher' wings.")
          )
        ),
        shiny::tabPanel(
          "Volcano (methods)",
          plotWidget(ns("volcanoMethods")),
          shiny::br(),
          tags$div(
            HTML("<b>Volcano plot for all statistical methods.</b> Simultaneous visualisation of volcano plots
                        of genes by multiple differential expression methods for the selected contrast.
                        Methods showing better statistical significance will show volcano plots with 'higher' wings.")
          )
        )
      )
    ),
    div(
      style = "max-height: 50vh",
      shiny::tabsetPanel(
        id = ns("tabs2"),
        shiny::tabPanel(
          "Table",
          tags$div(
            HTML("<b>Differential Expression Analysis.</b> Compare expression between
                        two conditions. Determine which genes are significantly downregulated or overexpressed in one of the groups.")
          ),
          shiny::br(),
          div(
            class = "row",
            div(
              class = "col-md-8",
              tableWidget(ns("genetable"))
            ),
            div(
              class = "col-md-4",
              tableWidget(ns("gsettable"))
            )
          )
        ),
        shiny::tabPanel(
          "Foldchange (all)",
          tableWidget(ns("fctable"))
        ),
        shiny::tabPanel(
          "FDR table",
          tableWidget(ns("FDRtable"))
        )
      )
    )
  )
}
