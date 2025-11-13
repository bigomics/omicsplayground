##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

upload_module_initial_settings_ui <- function(id) {
  ns <- shiny::NS(id)
  ## shiny::uiOutput(ns("UI"), fill = TRUE)

  div(
    style = "overflow: auto;",
    bslib::as_fill_carrier(),
    bslib::layout_columns(
      fill = FALSE,
      div(
        style = "display: flex; flex-direction: column; align-items: center; gap: 20px;",
        div(
          style = "margin-top: 0px; width: 40%;",
          bs_alert("This dialog will guide you to upload your own data into Omics Playground. You should prepare at the minimum two CSV files: an <b>expression.csv</b> file (containing your transcriptomics data) and a <b>samples.csv</b> file (containing your sample information). A third <b>contrasts.csv</b> file (describing your comparisons) is optional. Read more about data preparation <a href='https://omicsplayground.readthedocs.io/en/latest/dataprep/dataprep/'><u>here</u></a>.", closable = FALSE, translate = TRUE, html = TRUE)
        ),
        br(),
        div(
          p("Organism:", style = "text-align: left;   margin: 0 0 2px 0; font-weight: bold;"),
          shiny::selectInput(
            inputId = ns("selected_organism"),
            label = NULL,
            choices = playbase::SPECIES_TABLE$species_name,
            ## selected = 1,
            multiple = TRUE
          )
        ),
        div(
          p("Data type:", style = "text-align: left;   margin: 0 0 2px 0; font-weight: bold;"),
          shiny::selectInput(
            ns("selected_datatype"), NULL,
            choices = c(
              "RNA-seq",
              "scRNA-seq",
              "mRNA microarray",
              "proteomics",
              "other"
            )
          )
        )
        ## shiny::uiOutput(ns("upload_annot_ui"))
      )
    )
  )
}

upload_module_initial_settings_server <- function(
  id,
  upload_datatype,
  upload_organism,
  auth,
  new_upload
) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      observeEvent(new_upload(), {
        species <- playbase::SPECIES_TABLE$species_name
        if (!auth$options$ENABLE_ANNOT) {
          species <- setdiff(species, "No organism")
        }
        updateSelectizeInput(session, "selected_organism", choices = species, server = TRUE)
      })

      output$upload_annot <- renderUI({
        upload_annot_ui <- NULL
        if (auth$options$ENABLE_ANNOT) {
          upload_annot_ui <- fileInput2(
            ns("upload_annot_table"),
            shiny::tags$h4("Probe annotation (alpha):"),
            multiple = FALSE,
            accept = c(".csv")
          )
        }
        upload_annot_ui
      })

      # change upload_datatype to selected_datatype

      observeEvent(input$selected_datatype, {
        upload_datatype(input$selected_datatype)
      })

      # change upload_organism to selected_organism
      observeEvent(input$selected_organism, {
        upload_organism(input$selected_organism)
      })
    } ## end-of-server
  )
}
