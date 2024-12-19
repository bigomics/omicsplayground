##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

downloadButton2 <- function(outputId, label = "Download", class = NULL, ...) {
  aTag <- shiny::tags$a(
    id = outputId,
    class = paste("btn btn-default shiny-download-link", class),
    href = "", target = "_blank", download = NA,
    shiny::icon("file-csv"), label, ...
  )
}

LoadingInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings()
}

LoadingUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  board_header <- fillRow(
    flex = c(NA, NA, 1),
    shiny::div(
      id = "navheader-current-section",
      HTML("Load from Library &nbsp;"),
      shiny::actionLink(
        ns("module_info"), "",
        icon = shiny::icon("youtube"),
        style = "color: #ccc;"
      )
    ),
    shiny::div(shiny::uiOutput(ns("pgx_stats_ui")), id = "navheader-dataset-stats")
  )

  user_tabpanel <- shiny::tabPanel(
    "My Datasets",
    bslib::layout_columns(
      col_widths = 12,
      height = "calc(100vh - 180px)",
      uiOutput(ns("sharing_alert")),
      div(
        shiny::actionButton(
          ns("loadbutton"),
          label = "Load selected",
          icon = icon("file-import"),
          class = "btn btn-primary",
          width = NULL
        ),
        DatasetReportUI(id = ns("generate_report"))
      ),
      bslib::layout_columns(
        col_widths = c(8, 4),
        loading_table_datasets_ui(
          ns("pgxtable"),
          title = "Available datasets",
          info.text = "This table contains information about all available datasets within the platform. For each dataset, it reports a brief description as well as the total number of samples, genes, genesets, corresponding phenotypes and the creation date.",
          caption = "Table of datasets available in the platform.",
          height = c("calc(100vh - 340px)", 700),
          width = c("100%", "100%")
        ),
        loading_tsne_ui(
          ns("tsne"),
          title = "Signature t-SNE",
          info.text = "Scatter plot displaying the t-SNE clustering of the available contrasts. Each dot corresponds to a specific comparison.",
          info.methods = "t-SNE is a non-linear dimensionality reduction method that enables visualization of high-dimensional data in a low-dimensional space, typically 2D or 3D. Unlike linear dimensionality reduction techniques like PCA, t-SNE may separate data that is not linearly separable. Signatures/datasets that are clustered closer together, are more similar. Performed using the Rtsne R package [1].",
          info.references = list(
            list(
              "Krijthe JH (2015). Rtsne: T-Distributed Stochastic Neighbor Embedding using Barnes-Hut Implementation. R package version 0.17",
              "https://doi.org/10.32614/CRAN.package.Rtsne"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
          caption = "Similarity clustering of fold-change signatures colored by data sets using t-SNE.",
          height = c("calc(100vh - 340px)", "70vh"),
          width = c("auto", "100%")
        )
      ) ## end of 7fr-5fr
    )
  )

  public_tabpanel <- shiny::tabPanel(
    "Public Datasets",
    bslib::layout_columns(
      col_widths = 12,
      height = "calc(100vh - 180px)",
      bs_alert("This panel shows all <b>Public datasets</b>. You can select a public dataset and click <b>Import Dataset</b> to copy that dataset to your library for further analysis. The <b>Signature t-SNE</b> shows similarity clustering of fold-change signatures using t-SNE.", translate = FALSE, html = TRUE),
      bslib::layout_columns(
        col_widths = c(8, 4),
        loading_table_datasets_public_ui(
          ns("pgxtable_public"),
          title = "Public datasets",
          info.text = "This table shows available public datasets within the platform. For each dataset, it reports a brief description as well as the total number of samples, genes, gene sets (or pathways), corresponding phenotypes and the creation date.",
          caption = "Table with public datasets available in the platform.",
          ##height = c("calc(100vh - 330px)", 700),
          height = c("100%", 700),
          width = c("100%", "100%")
        ),
        loading_tsne_ui(
          ns("tsne_public"),
          title = "Signature t-SNE",
          info.text = "Scatter plot displaying the t-SNE clustering of the available contrasts. Each dot corresponds to a specific comparison.",
          info.methods = "t-SNE is a non-linear dimensionality reduction method that enables visualization of high-dimensional data in a low-dimensional space, typically 2D or 3D. Unlike linear dimensionality reduction techniques like PCA, t-SNE may separate data that is not linearly separable. Signatures/datasets that are clustered closer together, are more similar. Performed using the Rtsne R package [1].",
          info.references = list(
            list(
              "Krijthe JH (2015). Rtsne: T-Distributed Stochastic Neighbor Embedding using Barnes-Hut Implementation. R package version 0.17",
              "https://doi.org/10.32614/CRAN.package.Rtsne"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
          caption = "Similarity clustering of fold-change signatures colored by data sets using t-SNE.",
          ##height = c("calc(100vh - 330px)", 700),
          height = c("100%", 700),
          width = c("auto", "100%")
        )
      ) ## end of 7fr-5fr
    )
  ) ## end of Public tabPanel

  ## ------------------------------------------------------------------------

  ## disable/hide public tabpanel if public folder does not exists
  public_dir <- file.path(OPG, "data_public")
  if (!dir.exists(public_dir)) {
    public_tabpanel <- NULL
  }

  ## ============================ Board object ===========================
  div(
    class = "p-0",
    board_header,
    shiny::tabsetPanel(
      id = ns("tabs"),
      user_tabpanel,
      public_tabpanel
    )
  )
}

## ====================================================================
## ====================================================================
## ====================================================================


SharedDatasetsUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  tab_content <- bslib::layout_columns(
    col_widths = 12,
    height = "calc(100vh - 180px)",
    bs_alert("This Sharing panel shows <strong>received datasets</strong> that are not yet imported to your library, and your <strong>shared datasets</strong> that are still waiting to be accepted by the receiver. Please accept or refust each received file, and/or resend a message or cancel your shared datasets."),
    bslib::layout_columns(
      col_widths = 12,
      height = "calc(100vh - 180px)",
      uiOutput(ns("sharing_panel_ui"))
      ##      sharing_tabpanel
    )
  )

  div(
    class = "row",
    boardHeader(title = "Shared datasets", info_link = ns("loading_sharing")),
    ##    shiny::tabsetPanel(
    ##      id = ns("tabs1"),
    ##      shiny::tabPanel(
    ##        "Sharing",
    tab_content
    ##      ) ## tabPanel
    ##    ) ## tabsetPanel
  ) ## div
}
