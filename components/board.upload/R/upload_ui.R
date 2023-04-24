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
      bs_alert("In this panel, you can upload your data to the platform. The platform requires 3 data files as explained below: a data file containing counts/expression (counts.csv), a sample information file (samples.csv) and a file specifying the statistical comparisons as contrasts (contrasts.csv). NB Users can now create contrasts from the platform itself, so the contrasts.csv file is optional."),
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
      br(),
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
      bs_alert("Here, you can interactively create your comparisons (or so-called 'contrasts'). Choose a phenotype on the left, then create groups by dragging the conditions to the boxes of 'main' or 'control' group. Then click 'add comparison'. "),
      MakeContrastUI(ns("makecontrast"))
    ),
    shiny::tabPanel(
      "Compute",
      bs_alert("OK. We now have everything to compute your data. Please name your dataset and give a short description of the experiment. You can select/deselect some computation options but if you do not understand, it is safer to leave the defaults. If you are ready, hit 'Compute'. Computation will take 10-40 minutes depending on the size of your data and number of comparisons."),
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
