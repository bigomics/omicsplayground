##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

UploadInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    withTooltip(
      shiny::selectInput(ns("fa_contrast"), "Comparison:",
        choices = NULL
      ),
      "Select the comparison of interest.",
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
      bs_alert("In this panel, you can upload your data to the platform. The platform requires 3 data files as explained below: a data file containing counts/expression (counts.csv), a sample information file (samples.csv) and a file specifying the statistical comparisons (comparisons.csv). NB Users can now create comparisons from the platform itself, so the comparisons.csv file is optional."),
      br(),
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
          upload_plot_countstats_ui(
              id = ns("countStats"),
              title = "Count Stats",
              info.text = "Information about the uploaded counts.",
              caption = "Information about the uploaded counts.",
              height = c("75%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
          )
        ),
        div(
          class = "col-md-4",
          upload_plot_phenostats_ui(
              id = ns("phenoStats"),
              title = "Pheno Stats",
              info.text = "Information about the uploaded samples",
              caption = "Information about the uploaded samples.",
              height = c("75%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
          )
        ),
        div(
          class = "col-md-4",
          upload_plot_contraststats_ui(
              id = ns("contrastStats"),
              title = "Comparison Stats",
              info.text = "Information about the uploaded comparisons",
              caption = "Information about the uploaded comparisons.",
              height = c("75%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
          )
        )
      )
    ),
    shiny::tabPanel(
      "BatchCorrect",
      bs_alert("Omics data often suffers from batch effect due to experiments done on different days, using different machines or done at different institutes. This will often cause so-called batch effects. Batch correction can clean your data from these 'unwanted variation'. But be careful, batch correction can also be dangerous if not used carefully and can remove valuable real signal. Only adviced for advanced users!"),
      br(),
      shiny::fillCol(
        height = height,
        upload_module_batchcorrect_ui(ns("batchcorrect"))
      )
    ),
    shiny::tabPanel(
      "Comparisons",
      bslib::layout_column_wrap(
        width = 1,
        height = "calc(100vh - 200px)",
        heights_equal = "row",
        bs_alert(HTML("Here, you can interactively <b>create comparisons</b> (also called 'contrasts', 'groups'...). Choose a phenotype, then create groups by dragging conditions to the boxes of 'main' or 'control' group. Give the contrast a name (please keep it short!) and then click 'add comparison'. If you are feeling lucky, you can also try 'auto-comparisons'.")),
        upload_module_makecontrast_ui(ns("makecontrast"))
      )
    ),
    shiny::tabPanel(
      "Compute",
      bs_alert("OK. We now have everything to compute your data. Please name your dataset and give a short description of the experiment. You can select/deselect some computation options but if you do not understand, it is safer to leave the defaults. If you are ready, hit 'Compute'. Computation can take 10-40 minutes depending on the size of your data and number of comparisons."),
      br(),
      shiny::fillCol(
        height = height, ## width = 1200,
        upload_module_computepgx_ui(ns("compute"))
      )
    )
  )

  page_ui <- div(
    boardHeader(title = "Upload data", info_link = ns("module_info")),
    tabs
  )
  return(page_ui)
}
