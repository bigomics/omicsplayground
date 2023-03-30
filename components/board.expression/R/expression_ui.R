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
  imgH <- "35vh" ## height of images

  div(
    boardHeader(title = "Differential expression", info_link = ns("gx_info")),
    div(
      tagList(
        div(
          # style = "max-height:50vh;",
          shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel(
              "Plot",
              div(
                class = "row",
                div(
                  class = "col-md-3",
                  expression_plot_volcano_ui(ns("plots_volcano"),
                    label = "a",
                    height = c(imgH, TABLE_HEIGHT_MODAL),
                    width = c("auto", "100%")
                  ),
                ),
                div(
                  class = "col-md-3",
                  expression_plot_maplot_ui(
                    id = ns("plots_maplot"),
                    label = "b",
                    height = c(imgH, TABLE_HEIGHT_MODAL),
                    width = c("auto", "100%")
                  ),
                ),
                div(
                  class = "col-md-3",
                  expression_plot_barplot_ui(
                    id = ns("plots_barplot"),
                    label = "c",
                    height = c(imgH, TABLE_HEIGHT_MODAL),
                    width = c("auto", "100%")
                  ),
                ),
                div(
                  class = "col-md-3",
                  expression_plot_topfoldchange_ui(
                    id = ns("plots_topfoldchange"),
                    label = "d",
                    height = c(imgH, TABLE_HEIGHT_MODAL),
                    width = c("auto", "100%")
                  ),
                )
              ),
              tags$div(
                class = "caption",
                HTML("<b>Expression plots</b> associated with the selected contrast. <b>(a)</b> Volcano-plot plotting fold-change versuson significance the x and y axes, respectively. <b>(b)</b> MA-plot plotting signal intensity versus fold-change on the x and y axes, respectively. <b>(c)</b> Sorted barplot of the top diffentially expressed genes with largest (absolute) fold-change for selected contrast. <b>(d)</b> Sorted barplot of the differential expression of the selected gene across all contrasts.")
              )
            ),
            shiny::tabPanel(
              "Top genes",
              bslib::layout_column_wrap(
                width = 1,
                expression_plot_topgenes_ui(
                    id = ns("topgenes"),
                    label = "a",
                    height = c(imgH, TABLE_HEIGHT_MODAL),
                    width = c("auto", "100%")
                )
              )
#              tags$div(
#                class = "caption",
#                HTML("<b>Top differentially expressed genes.</b> Expression barplots of the top
#                      most differentially (both positively and negatively) expressed genes for
#                      the selected contrast.")
#              )
            ),
            shiny::tabPanel(
              "Volcano (all)",
              bslib::layout_column_wrap(
                width = 1,
                expression_plot_volcanoAll_ui(ns("volcanoAll"),
                  label = "a",
                  height = c(imgH, TABLE_HEIGHT_MODAL),
                  width = c("auto", "100%")
                )
              )
            ),
            shiny::tabPanel(
              "Volcano (methods)",
              bslib::layout_column_wrap(
                width = 1,
                expression_plot_volcanoMethods_ui(
                  id = ns("volcanoMethods"),
                  label = "a",
                  height = c(imgH, TABLE_HEIGHT_MODAL),
                  width = c("auto", "100%")
                )
              )
            )
          )
        ),
        shiny::br(),
        ##        shiny::br(),
        div(
          # style = "max-height: 50vh",
          shiny::tabsetPanel(
            id = ns("tabs2"),
            shiny::tabPanel(
              "Table",
              tags$div(
                class = "caption",                
                HTML("<b>Differential Expression Analysis.</b> Compare expression between two conditions. Determine which genes are significantly downregulated or overexpressed in one of the groups.")
              ),
              div(
                class = "row",
                div(
                  class = "col-md-8",
                  expression_table_genetable_ui(
                    ns("genetable"),
                    width = c("100%", "100%"),
                    height = c("400px", TABLE_HEIGHT_MODAL)
                  )
                ),
                div(
                  class = "col-md-4",
                  expression_table_gsettable_ui(
                    ns("gsettable"),
                    width = c("100%", "100%"),
                    height = c("400px", TABLE_HEIGHT_MODAL)
                  )
                )
              )
            ),
            shiny::tabPanel(
              "Foldchange (all)",
              tags$div(
                class = "caption",                
                HTML("<b>Differential expression (fold-change) across all contrasts.</b> The
                column `rms.FC` corresponds to the root-mean-square fold-change across all contrasts.")
              ),
              expression_table_fctable_ui(
                ns("fctable"),
                width = c("100%", "100%"),
                height = c("400px", TABLE_HEIGHT_MODAL)
              )
            ),
            shiny::tabPanel(
              "FDR table",
              expression_table_FDRtable_ui(
                ns("FDRtable"),
                width = c("100%", "100%"),
                height = c("400px", TABLE_HEIGHT_MODAL)
              )
            )
          )
        )
      )
    )
  )
}
