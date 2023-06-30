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
      HTML("Load dataset &nbsp;"),
      shiny::actionLink(
        ns("module_info"), "",
        icon = shiny::icon("youtube"),
        style = "color: #ccc;"
      )
    ),
    shiny::div(shiny::uiOutput(ns("pgx_stats_ui")), id = "navheader-dataset-stats")
  )

  user_tabpanel <- shiny::tabPanel(
    'User',
    bslib::layout_column_wrap(
      width = 1,
      heights_equal = "row",
      height = "calc(100vh - 180px)",
      uiOutput(ns("receive_pgx_alert")),
      bslib::layout_column_wrap(
        width = 1,
        style = htmltools::css(grid_template_columns = "7fr 5fr"),
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
          info.text = "Each dot corresponds to a specific comparison. Signatures/datasets that are clustered closer together, are more similar.",
          caption = "Similarity clustering of fold-change signatures colored by data sets using t-SNE.",
          height = c("calc(100vh - 340px)", "70vh"),
          width = c("auto",  "100%")
        )
      ), ## end of 7fr-5fr
      div(
        id = "load-action-buttons",
        # this button is needed to trigger download but should be hidden
        shiny::downloadLink(
          ns("download_pgx_btn"),
          label = "",
          icon = NULL,
          width = '0%'
        ),
        # this button is needed to trigger download but should be hidden
        shiny::downloadLink(
          ns("download_zip_btn"),
          label = "",
          icon = NULL,
          width = '0%'
        ),
        shiny::actionButton(
          ns("loadbutton"),
              label = "Load Dataset", icon = icon("file-import"),
          class = "btn btn-outline-primary"
        )
      ) ## end of buttons div
    )
  )

  public_tabpanel <- shiny::tabPanel(
    'Public',
    bslib::layout_column_wrap(
      width = 1,
      heights_equal = "row",
      height = "calc(100vh - 180px)",
      bs_alert(HTML("This panel shows all <b>Public datasets</b>. You can select a public dataset and click <b>Import Dataset</b> to copy that dataset to your library for further analysis. The <b>Signature t-SNE</b> shows similarity clustering of fold-change signatures using t-SNE.")),
      bslib::layout_column_wrap(
        width = 1,
        style = htmltools::css(grid_template_columns = "7fr 5fr"),
        loading_table_datasets_public_ui(
          ns("pgxtable_public"),
          title = "Public datasets",
          info.text = "This table shows available public datasets within the platform. For each dataset, it reports a brief description as well as the total number of samples, genes, gene sets (or pathways), corresponding phenotypes and the creation date.",
          caption = "Table with public datasets available in the platform.",
          height = c("calc(100vh - 330px)", 700),
          width = c("100%", "100%")
        ),
        loading_tsne_ui(
          ns("tsne_public"),
          title = "Signature t-SNE",
          info.text = "Each dot corresponds to a specific comparison/signature. Signatures that are clustered closer together, are more similar.",
          caption = "Similarity clustering of fold-change signatures colored by data sets using t-SNE.",
          height = c("calc(100vh - 330px)", 700),
          width = c("auto",  "100%")
        )
      ), ## end of 7fr-5fr
      div(
        id = "load-action-buttons",
        shiny::actionButton(
          ns("importbutton"),
          label = "Import dataset", icon = icon("file-import"),
              class = "btn btn-outline-primary"
        )
          ) ## end of buttons div
    ) ## end first layout_column_wrap
  ) ## end of Public tabPanel

  ## disable/hide public tabpanel if public folder does not exists
  public_dir <- file.path(OPG,"data_public")
  if(!dir.exists(public_dir)) {
    public_tabpanel <- NULL
  }

  ## return object
  div(
    class = "p-0",
    board_header,
    shiny::tabsetPanel(
      id = ns('tabs'),
      user_tabpanel,
      public_tabpanel
    )
  )

}
