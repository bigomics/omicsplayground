##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

UploadInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    withTooltip(shiny::selectInput(ns("fa_contrast"), "Contrast:",
                                   choices = NULL),
                "Select the contrast corresponding to the comparison of interest.",
                placement = "top"
    ),
    withTooltip(shiny::actionLink(ns("fa_options"), "Options",
                                  icon = icon("cog", lib = "glyphicon")),
                "Show/hide advanced options",
                placement = "top"
    ),
    shiny::br(),
    shiny::conditionalPanel(
      "input.fa_options % 2 == 1",
      ns = ns,
      shiny::tagList(
        withTooltip(
          shiny::checkboxInput(ns("fa_filtertable"),
                               "filter signficant (tables)",
                               FALSE),
          "Click to filter the significant entries in the tables."
        )
      )
    )
  )
}

UploadUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  tabs <- shiny::tabsetPanel(
    id = ns("tabs"),
    shiny::tabPanel(
      "Upload",
      div(
        class = "row",
        div(
          class = "col-md-3",
          shiny::sidebarPanel(
            width = "100%",
            fileInput2(ns("upload_files"),
                       shiny::h4("Choose files"),
                       multiple = TRUE, accept = c(".csv", ".pgx")
            ),
            shinyWidgets::prettySwitch(ns("load_example"), "Load example data"),
            shinyWidgets::prettySwitch(ns("advanced_mode"), "Batch correction (beta)")
          )
        ),
        div(
          class = "col-md-9",
          shiny::HTML(
            "<h4>User file upload</h4><p>Please prepare the data files
            in CSV format as listed below. It is important to name the files
            exactly as shown. The file format must be comma-separated-values
            (CSV) text. Be sure the dimensions, rownames and column names match
            for all files. You can download a zip file with example files here:
            EXAMPLEZIP. You can upload a maximum of <u>LIMITS</u>."
          )
        )
      ),
      div(
        class = "row",
        div(
          class = "col-md-4",
          shiny::plotOutput(ns("countStats")) %>% shinycssloaders::withSpinner()
        ),
        div(
          class = "col-md-4",
          shiny::plotOutput(ns("phenoStats")) %>% shinycssloaders::withSpinner()
        ),
        div(
          class = "col-md-4",
          shiny::plotOutput(ns("contrastStats")) %>% shinycssloaders::withSpinner()
        )
      )
    ),
    shiny::tabPanel(
      "BatchCorrect",
      shiny::fillCol(
        height = height,
        BatchCorrectUI(ns("batchcorrect"))
      )
    ),
    shiny::tabPanel(
      "Contrasts",
      shiny::fillCol(
                    height = 750,
                    flex = c(1,NA,NA,1),
                    shiny::fillRow(
                        flex = c(3,0.06,1.0),
                        shiny::fillCol(
                            flex = c(NA,NA,1.0),
                            shiny::h4("Create comparisons"),
                            ##p(help_text),
                            shiny::fillRow(
                                flex = c(1,4),
                                shiny::fillCol(
                                    flex = c(NA,NA,NA,NA,1),
                                    tipifyL(
                                        shiny::selectInput(ns("param"), "Phenotype:",
                                                    choices = NULL,
                                                    multiple = TRUE),
                                       "Select phenotype(s) to create conditions for your groups. Select <gene> if you want to split by high/low expression of some gene. Select <samples> if you want to group manually on sample names. You can select multiple phenotypes to create combinations."
                                    ),
                                    shiny::conditionalPanel(
                                        "input.param == '<gene>'", ns=ns,
                                        ##tipifyL(
                                        shiny::selectizeInput(ns("gene"), "Gene:", choices=NULL,
                                                       multiple=FALSE),
                                        ##"Select gene to divide your samples into high and low expression of that gene.")
                                    ),
                                    shiny::br(),
                                    tipifyL(
                                        shiny::textInput(ns("newname"), "Comparison name:",
                                                  placeholder="e.g. MAIN_vs_CONTROL"),
                                        "Give a name for your contrast as MAIN_vs_CONTROL, with the name of the main group first. You must keep _vs_ in the name to separate the names of the two groups."),
                                    shiny::br(),
                                    ## tipifyL(
                                    shiny::actionButton(ns("addcontrast"),
                                    "add comparison",
                                    icon=icon("plus"),
                                    class = "btn-outline-primary"),
                                    ##"After creating the groups, press this button to add the comparison to the table."a),
                                    shiny::br()
                                ),
                                withTooltip(
                                  shiny::uiOutput(ns("createcomparison"),
                                           style="font-size:13px; height: 280px; overflow-y: scroll;"),
                                  "Create comparisons by dragging conditions into the main or control groups on the right. Then press add comparison to add the contrast to the table.",
                                  placement="top", options = list(container = "body"))
                            )
                        ),
                        shiny::br(),
                        ##plotOutput(ns("pcaplot"), height="330px")
                        upload_plot_pcaplot_ui(
                          ns("pcaplot"),
                          height = c(320,700),
                          width = c("auto",800)
                        )

                        # plotWidget(ns("pcaplot"))
                    ),
                    shiny::h4("Contrast table"),
                    shiny::fillRow(
                        height = 24,
                        flex = c(NA,0.05,NA,NA,1),
                        withTooltip(
                            shiny::actionButton(ns("autocontrast"),
                            "add auto-contrasts",
                            icon=icon("plus"),
                            class="small-button btn-outline-primary"),
                            "If you are feeling lucky, try this to automatically create contrasts.",
                            placement="top", options = list(container = "body")
                        ),
                        shiny::br(),
                        shiny::div( shiny::HTML("<b>Strata:</b>"), style="padding: 4px 4px;"),
                        shiny::selectInput(ns("strata"), NULL, choices=NULL, width="120px"),
                        shiny::br()
                    ),
                    # shiny::br(),
                    ##shiny::tags$head(shiny::tags$style("table.dataTable.compact tbody th, table.dataTable.compact tbody td {padding: 0px 10px;}")),
                    ## this.style(ns("contrastTable"), "table.dataTable.compact tbody th, table.dataTable.compact tbody td {padding: 0px 10px;}"),
                    shiny::div(DT::dataTableOutput(ns("contrastTable")),
                        style="font-size:13px; height: 300px; margin-top: 20px;overflow-y: scroll;")
                )
      # shiny::uiOutput(ns("contrasts_UI"))
    ),
    shiny::tabPanel(
      "Compute",
      shiny::fillCol(
        height = height, ## width = 1200,
        ComputePgxUI(ns("compute"))
      )
    )
  )

  page_ui <- div(
    boardHeader(title = "Upload data", info_link = ns("module_info")),
    tabs
  )
  return(page_ui)
}
