##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

upload_module_initial_settings_ui <- function(id) {
  ns <- shiny::NS(id)

  div(
    style = "overflow: auto;",
    bslib::as_fill_carrier(),
    bslib::layout_columns(
      fill = FALSE,
      div(
        style = "display: flex; flex-direction: column; align-items: center; gap: 15px;",
        div(
          style = "margin-top: 0px; width: 45%;",
          bs_alert(HTML(tspan("This dialog will guide you to upload your own data into Omics Playground. You should prepare at the minimum two CSV files: an <b>expression.csv</b> file (containing your transcriptomics or proteomics data) and a <b>samples.csv</b> file (containing your sample information). A third <b>contrasts.csv</b> file (describing your comparisons) is optional. Read more about data preparation <a href='https://omicsplayground.readthedocs.io/en/latest/dataprep/dataprep/'><u>here</u></a>.")), closable = FALSE)
        ),
        br(),
        div(
          p("Organism:", style = "text-align: left;   margin: 0 0 2px 0; font-weight: bold;"),
          shiny::selectInput(
            inputId = ns("selected_organism"),
            label = NULL,
            choices = NULL,
            selected = 1,
            multiple = FALSE
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
        ),
        uiOutput(ns("upload_annot_table_ui"))
      )
    )
  )
}

upload_module_initial_settings_server <- function(
    id,
    upload_datatype,
    upload_organism,
    auth,
    new_upload) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      
      observeEvent( new_upload(), {
        # remove no organism
        species <- playbase::SPECIES_TABLE$species_name
        if (!auth$options$ENABLE_ANNOT) {
          species <- setdiff(species, "No organism")
        }
        shiny::updateSelectInput(
          session, 
          inputId = "selected_organism",
          choices = species
        )
      })
      
      output$upload_annot_table_ui <- shiny::renderUI({
          trigger <- new_upload()
          if (!auth$options$ENABLE_ANNOT) {
            return(NULL)
          }
          fileInput2(
            ns("upload_annot_table"),
            shiny::tags$h4("Probe annotation (alpha):"),
            multiple = FALSE,
            accept = c(".csv")
          )
      })

      # change upload_datatype to selected_datatype
      observeEvent(
        input$selected_datatype
      , {
        upload_datatype(input$selected_datatype)
        if (tolower(upload_datatype()) == "proteomics") {
          shiny.i18n::update_lang("proteomics", session)
        } else {
          shiny.i18n::update_lang("RNA-seq", session)
        }
      })

      # change upload_organism to selected_organism
      observeEvent(input$selected_organism, {
        upload_organism(input$selected_organism)
      })
      
    } ## end-of-server
  )
}
