##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

downloadButton2 <- function (outputId, label = "Download", class = NULL, ...) {
        aTag <- shiny::tags$a(id = outputId,
                       class = paste("btn btn-default shiny-download-link", class),
                       href = "", target = "_blank", download = NA, 
                       shiny::icon("file-csv"), label, ...)
}

LoadingInputs <- function(id) {
  ns <- shiny::NS(id)  ## namespace
  bigdash::tabSettings(
##    shiny::actionLink(ns("module_info"), "Tutorial", icon = shiny::icon("youtube")),
    shiny::hr(), shiny::br(),
    shiny::checkboxGroupInput(ns("flt_datatype"), "datatype", choices=""),
    shiny::checkboxGroupInput(ns("flt_organism"), "organism", choices="")
  )
}

LoadingUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace

    div(
      class = "p-1",
        uiOutput(ns("navheader")),
        br(), br(),

      ## table----------------            
      div(
          class = "row",
          div(
              class = "col-md-7",
              tableWidget(ns("pgxtable"))
          ),
          div(
              class = "col-md-5",
              loading_tsne_ui(ns("tsne"), height=c("65vh","70vh"))
          )
      ),
      br(),
          
      ## buttons----------------
      div( 
          id="load-action-buttons",
          shiny::actionButton(
              ns("deletebutton"), label="Delete dataset", icon=icon("trash"),
              class="btn btn-outline-danger-hover"
          ),
          shiny::downloadButton(
            ns("downloadpgx"), label="Download PGX", ##icon=icon("download"),
            class="btn btn-outline-dark-hover"
          ),
          downloadButton2(
            ns("downloadzip"), label="Download ZIP", icon=icon("file-archive"),
            class="btn btn-outline-dark-hover"
          ),
          shiny::actionButton(
            ns("loadbutton"), label="Load dataset", icon=icon("file-import"),
            class="btn btn-outline-primary"
          )          
      )
    )
}
