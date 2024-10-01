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
              "metabolomics (beta)" = "metabolomics"
              ## "scRNA-seq",
              ## "other"
            ),
            selected = DEFAULTS$datatype
          )
        ),
        shiny::uiOutput(ns("probe_type_ui")),
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
    useUploadWizard(ns),
    body
  )

  return(ui)
}


useUploadWizard <- function(ns) {

  counts_ui <- wizardR::wizard_step(
    step_title = tspan("Step 1: Upload counts", js = FALSE),
    step_id = "step_counts",
    upload_table_preview_counts_ui(
      ns("counts_preview")
    )
  )

  samples_ui <- wizardR::wizard_step(
    step_title = "Step 2: Upload samples",
    step_id = "step_samples",
    upload_table_preview_samples_ui(
      ns("samples_preview")
    )
  )

  contrasts_ui <- wizardR::wizard_step(
    step_title = "Step 3: Create comparisons",
    step_id = "step_comparisons",
    upload_table_preview_contrasts_ui(
      ns("contrasts_preview")
    )
  )

  ## batchcorrect_panel <- wizardR::wizard_step(
  ##   step_title = "BatchEffects",
  ##   step_id = "step_bc",
  ##   bslib::layout_columns(
  ##     col_widths = 12,
  ##     heights_equal = "row",
  ##     style = "margin-bottom: 20px",
  ##     upload_module_batchcorrect_ui(ns("batchcorrect")),
  ##   )
  ## )

  normalization_panel <- wizardR::wizard_step(
    step_title = "Step 4: QC/BC",
    step_id = "step_qc",
    upload_module_normalization_ui(ns("checkqc"))
  )

  compute_panel <- wizardR::wizard_step(
    step_title = "Compute!",
    step_id = "step_compute",
    upload_module_computepgx_ui(ns("compute"))
  )

  wizard <- wizardR::wizard(
    id = ns("upload_wizard"),
    width = 90,
    height = 75,
    modal = TRUE,
    style = "dots",
    lock_start = FALSE,
    ## initial_panel,
    counts_ui,
    samples_ui,
    contrasts_ui,
    normalization_panel,
    compute_panel,
    options = list(
      navigation = "buttons",
      finish = "Compute!"
    )
  )

  return(wizard)
}
