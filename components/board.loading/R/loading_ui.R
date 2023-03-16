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
    uiOutput(ns("navheader")),
    div(
      class = "row",
      div(
        class = "col-md-8",
        shiny::tabsetPanel(
          id = ns('tabs'),
          shiny::tabPanel(
            'User',
            loading_table_datasets_ui(
              ns("pgxtable"),
              height = c("65vh", 700),
              width = c("100%", "50%")
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
          shiny::tabPanel(
            'Shared',
            loading_table_datasets_shared_ui(
              ns("pgxtable_shared"),
              height = c("65vh", 700),
              width = c("100%", "50%")
            ),
            div(
              id = "load-action-buttons",
              shiny::actionButton(
                ns("importbutton"),
                label = "Import dataset", icon = icon("file-import"),
                class = "btn btn-outline-primary"
              )
            )
          )
        )
      ),
      div(
        class = "col-md-4",
        loading_tsne_ui(ns("tsne"),
        height = c("55vh", "60vh"),
        width = c("auto",  "100%")
        ) %>%
          tagAppendAttributes(
            style = 'padding-top: 61.5px;'
          )
      )
    )
  )
}
