##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

UploadUI <- function(id) {
  ns <- NS(id)

  species1 <- as.character(playbase::SPECIES_TABLE$species_name)
  species2 <- playbase::allSpecies.ORTHOGENE()
  AVAILABLE_SPECIES <- c(species1[1:4], intersect(species1, species2))

  body <- div(
    style = "overflow: auto;",
    bslib::as_fill_carrier(),
    bslib::layout_columns(
      fill = TRUE,
      div(
        style = "display: flex; flex-direction: column; align-items: center; gap: 20px; margin-bottom: 150px; margin-top: 120px;",
        div(
          style = "width: 40%;",
          bs_alert("To upload your own data, you should prepare at least two CSV files: an <b>expression.csv</b> file (containing your experiment data) and a <b>samples.csv</b> file (containing your sample information). A third <b>contrasts.csv</b> file (describing your comparisons) is optional. Read more about data preparation <a href='https://omicsplayground.readthedocs.io/en/latest/dataprep/dataprep/'><u>here</u></a>.", closable = FALSE, translate = TRUE, html = TRUE)
        ),
        br(),
        div(
          p("Data type:", style = "text-align: left; margin: 0 0 2px 0; font-weight: bold;"),
          shiny::selectInput(
            ns("selected_datatype"), NULL,
            choices = c(
              "RNA-seq",
              "proteomics",
              "mRNA microarray",
              "scRNA-seq",
              "other"
            )
          )
        ),
        div(
          p("Organism:", style = "text-align: left; margin: 0 0 2px 0; font-weight: bold;"),
          shiny::selectInput(
            inputId = ns("selected_organism"),
            label = NULL,
            choices = AVAILABLE_SPECIES,
            ## selected = 1,
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
    tagList(
      useUploadWizard(ns),
      body
    )
  )

  return(ui)
}


useUploadWizard <- function(ns) {
  ## initial_panel <- wizardR::wizard_step(
  ##   step_title = "Start",
  ##   step_id = "step_initial",
  ##   shiny::br(), shiny::br(),
  ##   ##        shinyWidgets::prettySwitch(ns("show_batchcorrection"), "Batch correction"),
  ##   ##        shinyWidgets::prettySwitch(ns("show_checkoutliers"), "Check outliers (beta)")
  ##   upload_module_initial_settings_ui(ns("initial"))
  ## )

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

  # comparisons_panel <- wizardR::wizard_step(
  #   step_title = "Comparison Builder",
  #   bslib::layout_columns(
  #     col_widths = 12,
  #     # height = "calc(100vh - 340px)",
  #     heights_equal = "row",
  #     upload_module_makecontrast_ui(ns("makecontrast")),
  #     bs_alert(HTML("Here, you can interactively <b>create comparisons</b> (also called 'contrasts'). Choose a phenotype, then create groups by dragging conditions to the boxes of the 'main' or 'control' group. Give the contrast a name (please keep it short!) and then click 'add comparison'. If you are feeling lucky, you can also try 'auto-comparisons'."))
  #   )
  # )

  batchcorrect_panel <- wizardR::wizard_step(
    step_title = "BatchEffects",
    step_id = "step_bc",
    bslib::layout_columns(
      col_widths = 12,
      heights_equal = "row",
      style = "margin-bottom: 20px",
      upload_module_batchcorrect_ui(ns("batchcorrect")),
      # bs_alert("Omics data often suffers from batch effect due to experiments done on different days, using different machines or done at different institutes. This will often cause so-called batch effects. Batch correction can clean your data from these 'unwanted variation'. But be careful, batch correction can also be dangerous if not used carefully and can remove valuable real signal. Only adviced for advanced users!")
    )
  )

  normalization_panel <- wizardR::wizard_step(
    step_title = "Step 4: QC/BC",
    step_id = "step_qc",
    upload_module_normalization_ui(ns("checkqc"))
  )

  compute_panel <- wizardR::wizard_step(
    step_title = "Compute!",
    step_id = "step_compute",
    # bs_alert("OK. We now have everything to compute your data. Please name your dataset and give a short description of the experiment. You can select/deselect some computation options but if you do not understand, it is safer to leave the defaults. If you are ready, hit 'Compute'. Computation can take 10-40 minutes depending on the size of your data and number of comparisons."),
    shiny::br(), shiny::br(),
    ## shinyWidgets::prettySwitch(ns("show_batchcorrection"), "Batch correction"),
    ## shinyWidgets::prettySwitch(ns("show_checkoutliers"), "Check outliers (beta)")
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

  wizard.save <- div(
    class = "p-0",
    ## div(
    ##   style = "position: fixed; right: 0px; width: 160px; margin-top: 10px;",
    ##   shinyWidgets::prettySwitch(ns("expert_mode"), "Expert mode")
    ## ),
    div(
      wizardR::wizard(
        id = ns("upload_wizard"),
        width = 90,
        height = 75,
        modal = TRUE,
        style = "dots",
        lock_start = TRUE,
        counts_ui,
        samples_ui,
        contrasts_ui,
        # comparisons_panel,
        # outliers_panel,
        # batchcorrect_panel,
        compute_panel,
        options = list(
          navigation = "buttons",
          finish = "Compute!"
        )
      )
    )
  )

  return(wizard)
}
