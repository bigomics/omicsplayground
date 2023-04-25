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
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    shiny::checkboxGroupInput(ns("flt_datatype"), "datatype", choices = ""),
    shiny::checkboxGroupInput(ns("flt_organism"), "organism", choices = "")
  )
}

LoadingUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  div(
    class = "p-0",
    fillRow(
        flex = c(NA, NA, 1),
        shiny::div(
          id = "navheader-current-section",
          HTML("Load dataset &nbsp;"),
          shiny::actionLink(
            ns("module_info"), "",
            icon = shiny::icon("info-circle"),
            style = "color: #ccc;"
          )
        ),
        shiny::div(shiny::uiOutput(ns("pgx_stats_ui")), id = "navheader-dataset-stats"),
    ),
    shiny::tabsetPanel(
      id = ns('tabs'),
      shiny::tabPanel(
        'User',
        bs_alert("This tab shows the available datasets within the platform. The table reports a brief description as well as the total number of samples, genes, gene sets (or pathways), corresponding phenotypes and the creation date. Select a dataset in the table and load the data by clicking the 'Load dataset' button."),
        div(
          class = "row",
          div(            
            class = "col-md-7",            
            loading_table_datasets_ui(
              ns("pgxtable"),
              title = "Available datasets",
              info.text = "This table contains information about all available datasets within the platform. For each dataset, it reports a brief description as well as the total number of samples, genes, genesets, corresponding phenotypes and the creation date.",
              caption = "Table of datasets available in the platform.",
              height = c("calc(100vh - 330px)", 700),
              width = c("100%", "100%")
            ),
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
            )
          ),
          div(
            class = "col-md-5",
            style = "height: 100%;",
            loading_tsne_ui(
              ns("tsne"),
              title = "Dataset explorer",
              info.text = "Each dot corresponds to a specific comparison. Signatures/datasets that are clustered closer together, are more similar.",
              caption = "Similarity clustering of fold-change signatures colored by data sets using t-SNE.",
              height = c("calc(100vh - 330px)", "70vh"),                           
              width = c("auto",  "100%")
            )
          )
        )
      ),

      shiny::tabPanel(
        'Shared',
        bs_alert("This tab shows all shared datasets. You can select a shared dataset and import that to your library for further analysis. You can also share any of your datasets to this shared folder from your library in the previous tab. Remember: sharing is caring!"),
        div(
          class = "row",
          div(
            class = "col-md-7",
            loading_table_datasets_shared_ui(
              ns("pgxtable_shared"),
              title = "Shared datasets",
              info.text = "This table shows available shared datasets within the platform. For each dataset, it reports a brief description as well as the total number of samples, genes, gene sets (or pathways), corresponding phenotypes and the creation date.",
              caption = "Table with shared datasets available in the platform.",
              height = c("calc(100vh - 330px)", 700),
              width = c("100%", "100%")
            ),
            div(
              id = "load-action-buttons",
              shiny::actionButton(
                ns("importbutton"),
                label = "Import dataset", icon = icon("file-import"),
                class = "btn btn-outline-primary"
              )
            )
          ),
          div(
            class = "col-md-5",
            loading_tsne_ui(
              ns("tsne_shared"),
              title = "Dataset explorer",
              info.text = "Each dot corresponds to a specific comparison/signature. Signatures that are clustered closer together, are more similar.",
              caption = "Similarity clustering of fold-change signatures colored by data sets using t-SNE.",
              height = c("calc(100vh - 330px)", 700),
              width = c("auto",  "100%")
            )
          )
        ) ## end of row
      ) ## end of tabPanel
    )
  )
}
