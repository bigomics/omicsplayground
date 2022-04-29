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

    shiny::tagList(
        shiny::fillCol(
          height = 750,
          flex = c(NA,NA,NA,1),
          ## shiny::fillRow(
          ##   height=115,
          ##   shiny::uiOutput(ns("valuebox1")),
          ##   shiny::uiOutput(ns("valuebox2")), 
          ##   shiny::uiOutput(ns("valuebox3"))
          ## ),
          uiOutput(ns("navheader")),
          tableWidget(ns("pgxtable")),

          div( id="load-action-buttons",
          shiny::fillRow(
            flex = c(NA,NA,NA,NA,1),
            withTooltip(
              shiny::actionButton(
                ns("loadbutton"), label="Load dataset", icon=icon("file-import"),
                class="btn btn-primary mx-2"),
              "Click to load the selected dataset.", placement="bottom"),
            withTooltip( shiny::downloadButton(
              ns("downloadpgx"), label="Download PGX", ##icon=icon("download"),
              class="btn btn-outline-primary mx-2")
             ,"Download PGX file (binary).", placement="bottom"),
            withTooltip( downloadButton2(
              ns("downloadzip"), label="Download ZIP", icon=icon("file-csv"),
              class="btn btn-outline-primary mx-2")
             ,"Download CSV files (counts.csv, samples.csv, contrasts.csv).",
              placement="bottom"),
            withTooltip( shiny::actionButton(
              ns("deletebutton"), label="Delete dataset", icon=icon("trash"),
              class="btn btn-outline-primary mx-2")
             ,"Delete the selected dataset.", placement="bottom"),
            br()
          )),
          br()
        )
    )
}
