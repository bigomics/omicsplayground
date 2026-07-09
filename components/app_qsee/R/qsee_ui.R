##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


qsee_ui <- function(id, height = "100%") {
  ns <- shiny::NS(id)

  bc_info <- HTML("<h4>Batch-effects analysis</h4>Batch correction can clean your data from 'unwanted variables'.\n")

  clust.infotext <-
    "Clustering of samples before ('uncorrected') and after the different batch-effect correction methods. After batch-effect correction clustering should improve. The silhouette score gives an indication of the clustering performance of the method."

  covariate.info <-
    "Analysis of variables by plotting their significance in correlation with the phenotype against their significance in correlation with a principal component (PC) vector. Strong model variables are situate 'top right'. Batch effect variables with high PC correlation but low phenotype correlation are on the 'top left'. A well-designed experiment shows strong model variables in PC1, else it may be a sign of significant batch-effects."

  pcc.info <- "PC analysis by covariate (class). The heights of the bars correspond to the relative contribution of that covariate to the three PC, as measured by an F-test. "

  clust.options <- tagList(
    shiny::radioButtons(ns("clust.plottype"), "Type;", c("tsne", "pca", "heatmap"))
  )

  covariate.options <- tagList(
    shiny::radioButtons(
      ns("covariate.plottype"), "Type;",
      c("Covariate plot", "Correlation heatmap")
    )
  )


  bslib::layout_columns(
    col_widths = c(2, 10),
    height = "calc(100vh - 200px)",
    heights_equal = "row",
    bslib::card(
      style = "background-color: #F7FAFD88;",
      bc_info,
      shiny::br(),
      shiny::br(),
      shiny::br(),
      withTooltip(
        shiny::actionLink(ns("adv_options"),
          "Advanced options",
          icon = icon("cog", lib = "glyphicon")
        ),
        "Toggle advanced options.",
        placement = "top"
      ),
      shiny::conditionalPanel(
        "input.adv_options % 2 == 1",
        ns = ns,
        shiny::selectizeInput(ns("tech_params"), "Technical covariates:",
          choices = "<none>", multiple = TRUE
        ),
        shiny::selectizeInput(ns("batch_params"), "Batch covariates:",
          choices = "<none>", multiple = TRUE
        ),
        shiny::actionButton(ns("recompute_button"), "Recompute",
          class = "btn-sm btn-primary mt-3", width = "100%"
        ),
      ),
      shiny::br(), shiny::br()
    ),
    bslib::layout_columns(
      width = 12,
      bslib::layout_columns(
        col_widths = 6,
        row_heights = c(3, 3),
        height = "calc(100vh - 200px)",
        heights_equal = "row",
        ##  shiny::plotOutput(ns("canvas"), width = "100%", height = height) %>% bigLoaders::useSpinner(),
        PlotModuleUI(
          ns("plot1"),
          title = "Clustering",
          info.text = clust.infotext,
          caption = clust.infotext,
          options = clust.options,
          height = c("100%", "70vh")
        ),
        PlotModuleUI(
          ns("plot2"),
          title = "PC components",
          info.text = pcc.info,
          caption = pcc.info,
          options = NULL,
          height = c("100%", "70vh")
        ),
        PlotModuleUI(
          ns("plot3"),
          title = "Covariate analysis",
          info.text = covariate.info,
          caption = covariate.info,
          options = covariate.options
        ),
        PlotModuleUI(
          ns("plot4"),
          title = "Statistics and score",
          #          info.text = info.text,
          #          caption = caption,
          options = NULL
        )
      )
    )
  )
}
