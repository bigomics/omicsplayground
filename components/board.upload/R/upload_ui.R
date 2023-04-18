##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

UploadInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    withTooltip(
      shiny::selectInput(ns("fa_contrast"), "Contrast:",
        choices = NULL
      ),
      "Select the contrast corresponding to the comparison of interest.",
      placement = "top"
    ),
    withTooltip(
      shiny::actionLink(ns("fa_options"), "Options",
        icon = icon("cog", lib = "glyphicon")
      ),
      "Show/hide advanced options",
      placement = "top"
    ),
    shiny::br(),
    shiny::conditionalPanel(
      "input.fa_options % 2 == 1",
      ns = ns,
      shiny::tagList(
        withTooltip(
          shiny::checkboxInput(
            ns("fa_filtertable"),
            "filter signficant (tables)",
            FALSE
          ),
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
          shiny::div(shiny::uiOutput(ns("upload_info")))
        )
      ),
      div(
        class = "row",
        div(
          class = "col-md-4",
          shiny::plotOutput(ns("countStats")) %>% bigLoaders::useSpinner()
        ),
        div(
          class = "col-md-4",
          shiny::plotOutput(ns("phenoStats")) # %>% bigLoaders::useSpinner()
        ),
        div(
          class = "col-md-4",
          shiny::plotOutput(ns("contrastStats")) # %>% bigLoaders::useSpinner()
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
      MakeContrastUI(ns("makecontrast"))
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
