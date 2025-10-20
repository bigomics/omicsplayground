##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

UploadUI <- function(id) {
  ns <- NS(id)

  body <- div(
    style = "overflow: auto;",
    bslib::as_fill_carrier(),
    bslib::layout_columns(
      fill = TRUE,
      div(
        style = "margin-top: 40px;",
        tabsetPanel(
          id = ns("upload_tabs"),
          type = "tabs",
          # Upload your own data tab
          tabPanel(
            "Upload your data",
            value = "upload",
            div(
              style = "display: flex; flex-direction: column; align-items: center; gap: 20px; margin-bottom: 150px; margin-top: 40px;",
              div(
                style = "width: 60%;",
                bs_alert("To upload your own data, you should prepare at least two CSV files: an <b>counts.csv</b> file (containing your experiment data) and a <b>samples.csv</b> file (containing your sample information). A third <b>contrasts.csv</b> file (describing your comparisons) is optional. Read more about data preparation <a href='https://omicsplayground.readthedocs.io/en/latest/dataprep/dataprep/'><u>here</u></a>.", closable = FALSE, translate = TRUE, html = TRUE)
              ),
              br(),
              div(
                p("Data type:", style = "text-align: left; margin: 0 0 2px 0; font-weight: bold;"),
                shiny::selectInput(
                  ns("selected_datatype"), NULL,
                  choices = c(
                    "RNA-seq",
                    "mRNA microarray",
                    "proteomics",
                    "scRNA-seq",
                    "metabolomics (beta)" = "metabolomics",
                    "multi-omics (beta)" = "multi-omics"
                  ),
                  selected = DEFAULTS$datatype,
                  width = "400px"
                ),
                shiny::uiOutput(ns("proteomics_subtype_ui"))
              ),
              div(
                p("Organism:", style = "text-align: left; margin: 0 0 2px 0; font-weight: bold;"),
                shiny::selectInput(
                  inputId = ns("selected_organism"),
                  label = NULL,
                  choices = NULL,
                  multiple = FALSE,
                  width = "400px"
                )
              ),
              br(),
              shiny::actionButton(
                ns("start_upload"),
                "Start upload",
                class = "btn-primary"
              ),
              br()
            )
          ),

          # Retrieve public datasets tab
          tabPanel(
            "Retrieve public datasets",
            value = "public",
            div(
              style = "display: flex; flex-direction: column; align-items: center; gap: 20px; margin-bottom: 150px; margin-top: 40px;",
              div(
                style = "width: 60%;",
                bs_alert("The 'Retrieve public data' functionality in OPG allows to query large repositories of publically available data to retrieve a specific dataset. You will need the dataset's unique identifier, such as the GEO ID. At the moment, <a href='https://www.ncbi.nlm.nih.gov/geo/' target='_blank'><u>GEO</u></a>, <a href='https://jhubiostatistics.shinyapps.io/recount3-study-explorer/' target='_blank'><u>recount3</u></a>, and <a href='https://www.ebi.ac.uk/biostudies/arrayexpress' target='_blank'><u>ArrayExpress</u></a> repositories are queried.", closable = FALSE, translate = TRUE, html = TRUE)
              ),
              br(),
              div(
                p("Data type:", style = "text-align: left; margin: 0 0 2px 0; font-weight: bold;"),
                shiny::selectInput(
                  ns("selected_datatype_public"), NULL,
                  choices = c(
                    "RNA-seq",
                    "mRNA microarray"
                  ),
                  selected = "RNA-seq",
                  width = "400px"
                )
              ),
              div(
                p("Organism:", style = "text-align: left; margin: 0 0 2px 0; font-weight: bold;"),
                shiny::selectInput(
                  inputId = ns("selected_organism_public"),
                  label = NULL,
                  choices = NULL,
                  multiple = FALSE,
                  width = "400px"
                )
              ),
              div(
                p("Dataset identifier:", style = "text-align: left; margin: 0 0 2px 0; font-weight: bold;"),
                shiny::textInput(
                  ns("dataset_identifier"),
                  label = NULL,
                  placeholder = "e.g., GSE12345, E-GEOD-12345",
                  width = "400px"
                ),
                div(
                  style = "font-size: 0.85em; color: #666; margin-top: 5px;",
                  "Supported: GEO (GSE), ArrayExpress (E-), recount3 datasets"
                )
              ),
              br(),
              shiny::actionButton(
                ns("start_search"),
                "Retrieve dataset",
                class = "btn-primary"
              ),
              br()
            )
          )
        )
      )
    )
  )

  ui <- div(
    boardHeader(title = "Upload New", info_link = ns("upload_info")),
    uiOutput(ns("upload_wizard")),
    body
  )
  return(ui)
}
