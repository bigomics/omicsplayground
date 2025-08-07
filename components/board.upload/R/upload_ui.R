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
        style = "display: flex; flex-direction: column; align-items: flex-start; gap: 20px; margin-bottom: 150px; margin-top: 250px;",

        div(
          style = "display: flex; flex-direction: row; align-items: center; gap: 20px;",
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] <= 0", ns("public_data_opts")),
            shiny::actionButton(ns("show_upload_opts"), "Upload your data", class = "btn-primary upload-btn-large")
          ),
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] <= 0", ns("show_upload_opts")),
            shiny::actionButton(ns("public_data_opts"), "Retrieve public data", class = "btn-primary upload-btn-large")
          )
        ),
        br(),

        tags$style(HTML("
          .upload-btn-large {
            font-size: 1.6em;
            padding: 20px 45px;
            min-width: 260px;
            border-radius: 8px;
            margin-left: 20vw;
            display: block;
            box-shadow: 0 2px 16px rgba(0,0,0,0.13);
          }
          .upload-btn-wrapper {
            width: 100%;
            display: flex;
            justify-content: flex-start;
            margin: 60px 5vw 0 0;
          }
        ")),

        shiny::conditionalPanel(
          condition = sprintf("input['%s'] > 0", ns("show_upload_opts")),
          div(
            style = "display: flex; flex-direction: column; align-items: center; gap: 20px; margin-bottom: 150px; margin-top: 120px;",
            div(
              style = "width: 60%;",
              bs_alert("To upload your own data, you should prepare at least two CSV files: an <b>counts.csv</b> file (containing your experiment data) and a <b>samples.csv</b> file (containing your sample information). A third <b>contrasts.csv</b> file (describing your comparisons) is optional. Read more about data preparation <a href='https://omicsplayground.readthedocs.io/en/latest/dataprep/dataprep/'><u>here</u></a>.",
                closable = FALSE, translate = TRUE, html = TRUE)
            ),
            br(),
            div(
              p("Data type:", style = "text-align: left; margin: 0 0 2px 0; font-weight: bold;"),
              shiny::selectInput(
                ns("selected_datatype"), NULL,
                choices = c("RNA-seq", "mRNA microarray", "proteomics", "scRNA-seq",
                  "metabolomics (beta)" = "metabolomics", "multi-omics (beta)" = "multi-omics"),
                selected = DEFAULTS$datatype,
                width = "400px"
              )
            ),
            div(
              p("Organism:", style = "text-align: left; margin: 0 0 2px 0; font-weight: bold;"),
              shiny::selectInput(
                inputId = ns("selected_organism"), label = NULL,
                choices = NULL, multiple = FALSE, width = "400px"
              )
            ),
            br(),
            shiny::actionButton(ns("start_upload"), "Start upload", class = "btn-primary"),
            br()
          )
        ),
        
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] > 0", ns("public_data_opts")),
          div(
            style = "display: flex; flex-direction: column; align-items: center; gap: 20px; margin-bottom: 150px; margin-top: 120px;",
            div(
              style = "width: 60%;",
              bs_alert("The 'Retrieve public data' functionality in OPG allows to query large repositories of publically available data to retrieve a specific dataset. You will need the dataset's unique identifier, such as the GEO ID. At the moment, GEO, ReCount, and ArrayExpress repositories are queried.",
                closable = FALSE, translate = TRUE, html = TRUE)
            ),
            br(),
            div(
              p("Dataset identifier:", style = "text-align: left; margin: 0 0 2px 0; font-weight: bold;"),
              shiny::textInput(ns("dataset_identifier"), label = NULL, width = "400px")
            ),
            br(),
            shiny::actionButton(ns("start_search"), "Start search & download", class = "btn-primary"),
            br()
          )
        ),

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
