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
    class = "p-1",
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
        div(
          class = "row",
          div(            
            class = "col-md-7",            
            loading_table_datasets_ui(
              ns("pgxtable"),
              title = "Available datasets",
              info.text = "This table contains information about all available datasets within the platform. For each dataset, it reports a brief description as well as the total number of samples, genes, genesets, corresponding phenotypes and the creation date.",
              caption = "Table of datasets available in the platform.",
              height = c("calc(100vh - (240px + 70px))", 700),
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
              height = c("calc(100vh - (240px + 70px))", "70vh"),                           
              width = c("auto",  "100%")
            )
          )
        )
      ),

      shiny::tabPanel(
        'Shared',
        div(
          class = "row",
          div(
            class = "col-md-7",
            loading_table_datasets_shared_ui(
              ns("pgxtable_shared"),
              title = "Shared datasets",
              info.text = "This table contains a general information about all available datasets within the platform. For each dataset, it reports a brief description as well as the total number of samples, genes, gene sets (or pathways), corresponding phenotypes and the creation date.",
              caption = "Table with public datasets available in the platform.",
              height = c("65vh", 700),
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
              info.text = "Each dot corresponds to a specific comparison. Signatures/datasets that are clustered closer together, are more similar.",
              caption = "Similarity clustering of fold-change signatures colored by data sets using t-SNE.",
              height = c("65vh", "70vh"),
              width = c("auto",  "100%")
            )
          )
        )
      )
    )
  )
}
