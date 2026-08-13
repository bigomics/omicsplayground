##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Was: board.upload/upload_module_batchcorrect.R


qsee_bsee_inputs <- function(id) {
  ns <- shiny::NS(id)

  bigdash::tabSettings(
    shiny::selectizeInput(ns("main_param"), "Main parameter:",
      choices = NULL, multiple = FALSE
    ),
    shiny::selectizeInput(ns("tech_params"), "Technical covariates:",
      choices = "<none>", multiple = TRUE
    ),
    shiny::selectizeInput(ns("batch_params"), "Batch covariates:",
      choices = "<none>", multiple = TRUE
    ),
    shiny::actionButton(ns("recompute_button"), "Recompute",
      class = "btn-sm btn-primary mt-3", width = "100%"
    )
  )
}

qsee_bsee_ui <- function(id) {
  ns <- shiny::NS(id)

  clust.infotext <-
    "Clustering of samples before ('uncorrected') and after the different batch-effect correction methods. After batch-effect correction clustering should improve. The silhouette score gives an indication of the clustering performance of the method."

  covariate.info <-
    "Analysis of variables by plotting their significance in correlation with the phenotype against their significance in correlation with a principal component (PC) vector. Strong model variables are situate 'top right'. Batch effect variables with high PC correlation but low phenotype correlation are on the 'top left'. A well-designed experiment shows strong model variables in PC1, else it may be a sign of significant batch-effects."

  pcc.info <- "PC analysis by covariate (class). The heights of the bars correspond to the relative contribution of that covariate to the three PC, as measured by an F-test. "

  clust.options <- tagList(
    shiny::selectInput(ns("clust.colorby"), "Color by:", choices = NULL, multiple = FALSE),
    shiny::checkboxInput(ns("show_labels"), "Show sample labels", value = FALSE)
  )

  batchcorrect_panel <-
    bslib::navset_tab(
      bslib::nav_panel(
        title = "Analysis",
        height = "100%",
        bslib::layout_columns(
          col_widths = 6,
          row_heights = c(1,1),
          height = "calc(100vh - 140px)",
          heights_equal = "row",
          PlotModuleUI(
            ns("plot1"),
            title = "Clustering",
            info.text = clust.infotext,
            caption = clust.infotext,
            plotlib = "plotly",
            options = clust.options,
            height = c("100%", "70vh")
          ),
          bslib::card(
            full_screen = TRUE,
            bslib::card_header("Heatmap"),
            qsee_plotly_hm_grid_ui(ns, "bsee", n = 6L, ncol = 3L)
          ),
          PlotModuleUI(
            ns("plot3"),
            title = "Statistics and score",
            plotlib = "plotly"
          ),
          PlotModuleUI(
            ns("plot4"),
            title = "Covariate correlation",
            info.text = covariate.info,
            caption = covariate.info,
            plotlib = "plotly"
          )
        )
      ),  ## end of Analysis panel (was Panel2)
      bslib::nav_panel(
        title = "PCVA",
        bslib::layout_columns(
          col_widths = 6,
          height = "calc(100vh - 140px)",
          bslib::layout_columns(
            col_widths = 12,
            row_heights = c(1,1),
            PlotModuleUI(
              ns("plot5"),
              title = "PVCA by phenotype",
              info.text = pcc.info,
              caption = pcc.info,
              plotlib = "plotly",
              height = c("100%", "70vh")
            ),
            PlotModuleUI(
              ns("plot6"),
              title = "PVCA by component",
              info.text = pcc.info,
              caption = pcc.info,
              plotlib = "plotly",
              height = c("100%", "70vh")
            )
          ),
          PlotModuleUI(
            ns("plot7"),
            title = "Covariate analysis",
            info.text = covariate.info,
            caption = covariate.info,
            plotlib = "plotly",
            height = c("100%", "70vh")
          )
        )
      )  ## end of Before vs After panel (was Panel1)
    )

  OmicsBoardUI(
    id = ns("board"),
    title = "Batch-effects",
    board_visibility_probe(ns),
    batchcorrect_panel
  )
}
