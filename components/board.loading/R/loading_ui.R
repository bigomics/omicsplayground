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
              width = c("100%", "50%")
            ),
            div(
              id = "load-action-buttons",
              shiny::actionButton(
                ns("deletebutton"),
                label = "Delete dataset", icon = icon("trash"),
                class = "btn btn-outline-danger"
              ),
              shiny::downloadButton(
                ns("downloadpgx"),
                label = "Download PGX",
                class = "btn btn-outline-dark"
              ),
              downloadButton2(
                ns("downloadzip"),
                label = "Download ZIP", icon = icon("file-archive"),
                class = "btn btn-outline-dark"
              ),
              shiny::actionButton(
                ns("loadbutton"),
                label = "Load dataset", icon = icon("file-import"),
                class = "btn btn-outline-primary"
              )
            )
          ),
          div(
            class = "col-md-5",
            loading_tsne_ui(ns("tsne"),
              height = c("65vh", "70vh"),
              width = c("auto",  "100%")
            ) %>%
              tagAppendAttributes(
                ##style = 'padding-top: 61.5px;'
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
          ),
          div(
            class = "col-md-5",
            loading_tsne_ui(ns("tsne_shared"),
              height = c("65vh", "70vh"),
              width = c("auto",  "100%")
            ) %>%
              tagAppendAttributes(
                ## style = 'padding-top: 61.5px;'
              )
          )          
        )
      )
    )
  )
}
