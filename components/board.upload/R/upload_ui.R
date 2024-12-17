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
        style = "display: flex; flex-direction: column; align-items: center; gap: 20px; margin-bottom: 150px; margin-top: 120px;",
        div(
          style = "width: 40%;",
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
              ## "scRNA-seq",
              ## "other"
            ),
            selected = DEFAULTS$datatype
          )
        ),
        ##        shiny::uiOutput(ns("probe_type_ui")),
        div(
          p("Organism:", style = "text-align: left; margin: 0 0 2px 0; font-weight: bold;"),
          shiny::selectInput(
            inputId = ns("selected_organism"),
            label = NULL,
            choices = NULL,
            multiple = FALSE
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
    )
  )

  ui <- div(
    boardHeader(title = "Upload New", info_link = ns("upload_info")),
    uiOutput(ns("upload_wizard")),
    body
  )
  return(ui)
}
