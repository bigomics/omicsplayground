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
    shiny::tagList(
        shinyBS::tipify( shiny::actionLink(ns("module_info"), "Tutorial", icon = shiny::icon("youtube")),
               "Show more information about this module.")
    )
}

LoadingUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace

    close_session <- shiny::span()
    if(getOption("OMICS_TEST", FALSE)){
        close_session <- shiny::actionButton(ns("close"), "close")
    }

    shiny::tagList(
        close_session,
        shiny::fillCol(
            height = 750,
            shiny::tabsetPanel(
                id = ns("tabs"),
                shiny::tabPanel("Datasets",
                    shiny::fillCol(
                        height = 750,
                        flex = c(NA,1),
                        shiny::fillRow(
                            height=115,
                            shiny::uiOutput(ns("valuebox1")),
                            shiny::uiOutput(ns("valuebox2")), 
                            shiny::uiOutput(ns("valuebox3"))
                        ),
                        shiny::fillRow(
                            flex = c(1,0.1,4.5),
                            shiny::wellPanel(
                                shiny::tagList(
                                    shinyalert::useShinyalert(),
                                    shiny::p(shiny::strong("Dataset info:")),
                                    shiny::div(shiny::htmlOutput(ns("dataset_info")), id="datainfo"),
                                    shiny::br(),
                                    shiny::conditionalPanel(
                                        "output.rowselected != 0", ns=ns,
                                        shinyBS::tipify( shiny::actionButton(ns("loadbutton"),label="Load",class="load-button"),
                                        "Click to load the selected dataset.", placement="bottom"),
                                        shinyBS::tipify( shiny::downloadButton(
                                            ns("downloadpgx"), label=NULL, ## icon=icon("download"),
                                            style='padding:2px 1px 1px 1px; font-size:140%; width:30px;'
                                        ),"Download PGX file (binary).", placement="bottom"),
                                        shinyBS::tipify( downloadButton2(
                                            ns("downloadzip"), label=NULL, icon=icon("file-csv"),
                                            style='padding:2px 1px 1px 1px; font-size:140%; width:30px;'
                                        ),"Download CSV files (counts.csv, samples.csv, contrasts.csv).",
                                        placement="bottom"),
                                        shinyBS::tipify( shiny::actionButton(
                                            ns("deletebutton"), label=NULL, icon=icon("trash"),
                                            style='padding:2px 1px 1px 1px; font-size:140%; color: #B22222; width:30px;'
                                        ),"Delete the selected dataset.", placement="bottom")
                                    ),
                                    shiny::br(),shiny::br(),
                                    shinyBS::tipify( shiny::actionLink(ns("showfilter"), "show filters",
                                                                    icon=icon("cog", lib = "glyphicon")),
                                        "Show dataset filters.", placement="top"),
                                    shiny::br(),br(),
                                    shiny::conditionalPanel(
                                        "input.showfilter % 2 == 1", ns=ns,
                                        shiny::uiOutput(ns("dataset_filter"))
                                    )
                                )
                            ),
                            shiny::br(), 
                            tableWidget(ns("pgxtable"))
                        )
                    )
                ),
                shiny::tabPanel("Upload data",
                    UploadModuleUI(ns("upload_panel"))
                )
            )
        )
    )
}
