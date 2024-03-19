##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


UploadUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  counts_ui <- wizardR::wizard_step(
    step_title = "Counts",
    upload_table_preview_counts_ui(
      ns("counts_preview")
      )
  )
  
  samples_ui <- wizardR::wizard_step(
    step_title = "Samples",
    upload_table_preview_samples_ui(
      ns("samples_preview")
    )
  )

  contrasts_ui <- wizardR::wizard_step(
    step_title = "Comparison",
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
    bslib::layout_columns(
      col_widths = 12,
      heights_equal = "row",
      upload_module_batchcorrect_ui(ns("batchcorrect")),
      #bs_alert("Omics data often suffers from batch effect due to experiments done on different days, using different machines or done at different institutes. This will often cause so-called batch effects. Batch correction can clean your data from these 'unwanted variation'. But be careful, batch correction can also be dangerous if not used carefully and can remove valuable real signal. Only adviced for advanced users!")      
    )
  )

  outliers_panel <- wizardR::wizard_step(
    step_title = "QC/BC",
    bslib::layout_columns(
      col_widths = 12,
      # height = "calc(100vh - 340px)",
      heights_equal = "row",
      upload_module_outliers_ui(ns("checkqc")),
      #bs_alert("Check for normalization, outliers and batch-effects.")
    )
  )
  
  compute_panel <- wizardR::wizard_step(
    step_title = "Dataset description",
    #bs_alert("OK. We now have everything to compute your data. Please name your dataset and give a short description of the experiment. You can select/deselect some computation options but if you do not understand, it is safer to leave the defaults. If you are ready, hit 'Compute'. Computation can take 10-40 minutes depending on the size of your data and number of comparisons."),
    shiny::br(), shiny::br(),
    
##        shinyWidgets::prettySwitch(ns("show_batchcorrection"), "Batch correction"),
##        shinyWidgets::prettySwitch(ns("show_checkoutliers"), "Check outliers (beta)")
      upload_module_computepgx_ui(ns("compute"))
  )

  
review_panel <- wizardR::wizard_step(
  step_title = "Review and compute",
  shiny::br(), shiny::br(),
  shiny::htmlOutput(ns("input_recap")),
  shiny::br(), shiny::br(),
  shiny::fluidRow(
    bslib::layout_columns(
      col_widths = c(4,4,4),
      upload_plot_countstats_ui(
        id = ns("countStats"),
        title = "Count Stats",
        info.text = "Information about the uploaded counts.",
        caption = "Information about the uploaded counts.",
        height = c("auto", "100%"),
        width = c("auto", "100%")
      ),
      upload_plot_phenostats_ui(
        id = ns("phenoStats"),
        title = "Pheno Stats",
        info.text = "Information about the uploaded samples",
        caption = "Information about the uploaded samples.",
        height = c("auto", "100%"),
        width = c("auto", "100%")
      ),
      upload_plot_contraststats_ui(
        id = ns("contrastStats"),
        title = "Comparison Stats",
        info.text = "Information about the uploaded comparisons",
        caption = "Information about the uploaded comparisons.",
        height = c("auto", "100%"),
        width = c("auto", "100%")
      )
    )
  )
)
    div(
    class = "p-0",
    # board_header,
    div(
      style = "position: fixed; right: 0px; width: 160px; margin-top: 10px;",
      shinyWidgets::prettySwitch(ns("expert_mode"), "Expert mode"),
    ),
    div(
    wizardR::wizard(
      id = ns("upload_wizard"),
      width = 90,
      height = 80,
      modal = TRUE,
      style = "progress",
      lock_start = TRUE,
      counts_ui,
      samples_ui,
      contrasts_ui,
      # comparisons_panel,
      # outliers_panel,
      # batchcorrect_panel,
      compute_panel,
      review_panel,
      options = list(
        navigation = "buttons",
        finish = "Compute!"
      )
    )
    )
  )
}
