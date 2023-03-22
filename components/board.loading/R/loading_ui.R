##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
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
        shiny::div(selector_default(ns("hide_caption"), label = "Show captions"))
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
              height = c("65vh", 700),
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
            loading_tsne_ui(ns("tsne"),
              height = c("65vh", "70vh"),
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
            loading_tsne_ui(ns("tsne_shared"),
              height = c("65vh", "70vh"),
              width = c("auto",  "100%")
            )
          )
        )
      )
    )
  )
}
