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
            "filter significant (tables)",
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

  board_header <- fillRow(
    flex = c(NA, NA, 1),
    shiny::div(
      id = "navheader-current-section",
      HTML("Create New Dataset &nbsp;"),
      shiny::actionLink(
        ns("module_info"), "",
        icon = shiny::icon("youtube"),
        style = "color: #ccc;"
      )
    )
  )

  upload_select_db <- shiny::tabPanel(
    "Select Organism",
    bslib::layout_column_wrap(
      width = 1,
      heights_equal = "row",
      height = "calc(100vh - 180px)",
      # add a drop down selector for organism
      shiny::div(
        style = "display: flex; justify-content: center; align-items: center; margin-top: 100px;",
        shiny::div(
          style = "text-align: center;",
          h3(shiny::HTML("<b>Select the organism:</b>")),
          div(
            style = "margin-top: 30px; padding-left: 70px; text-align: center;",
            shiny::selectInput(
              ns("selected_organism"),
              NULL,
              # restrict to ensembl species, as we are validating them in the first place
              choices = playbase::SPECIES_TABLE$species_name[which(playbase::SPECIES_TABLE$mart == "ensembl")],
              selected = NULL,
              multiple = FALSE
            )
          ),
          shiny::div(
            style = "margin-top: 20px;text-align: center;",
            shiny::actionButton(ns("proceed_to_upload"), "Next",
              icon = icon("arrow-right"),
              class = "btn btn-success"
            )
          ),
          div(
            style = "margin-top: 120px",
            h3("Need a dataset to try?")
          ),
          shiny::div(
            style = "margin-top: 30px",
            shiny::downloadButton(
              ns("downloadExampleData2"),
              width = "220px",
              icon = icon("download"),
              label = "Download example data",
              class = "btn-outline-primary",
              style = "margin-right: 10px;"
            ),
            shiny::actionButton(
              ns("load_example2"),
              width = "220px",
              icon = icon("table"),
              label = "Use example data",
              class = "btn-outline-primary",
              style = "margin-left: 10px;"
            )
          )
        )
      )
    )
  )


  upload_panel <- shiny::tabPanel(
    "Upload",
    bslib::layout_columns(
      col_widths = 12,
      height = "calc(100vh - 180px)",
      bslib::layout_columns(
        col_widths = c(4, 8),        
        bslib::card(
          height = "100%",
          width = "100%",
          style = "background-color: #f7fafd;",
          shiny::selectInput(
            ns("selected_organism"),
            h4("1. Select organism:", class='mb-0'),
            ## restrict to ensembl species, as we are validating them in the first place
            choices = playbase::SPECIES_TABLE$species_name[which(playbase::SPECIES_TABLE$mart == "ensembl")],
            selected = NULL,
            multiple = FALSE,
            width = "90%"            
          ),
##          shiny::br(), 
          div(
            style = 'display: none',
            fileInput2(
              ns("upload_files"),              
              shiny::h4("Choose files", class='mb-0'),
              multiple = TRUE,
              buttonClass = "btn-primary",
              accept = c(".csv", ".pgx")
            )
          ),
          shiny::h4("2. Choose data:", class='mt-2'),
          div(
            style = "margin-left: 0px; margin-right: auto;",
            shiny::actionButton(
              ns("upload_files_btn"),
              width = 'auto',
              icon = icon("table"),
              label = "Upload files",
              class = "btn-primary",
              style = "margin-left: 0px;"
            ),
            shiny::actionButton(
              ns("load_example"),
              width = 'auto',                
              icon = icon("table"),
              label = "Use example data",
              class = "btn-outline-info",
              style = "margin-left: 0px;"
            )
          )
        ),
        bslib::card(
          shiny::h4("How to upload your files."),
          shiny::uiOutput(ns("upload_info"), class='mt-3 mb-1'),
          shiny::div(
            ##style = "margin-top: 10px; position: absolute; bottom: 47px;",
            style = "margin-top: 25px;",            
            shiny::downloadButton(
              ns("downloadExampleData"),
              width = "200px",
              icon = icon("download"),
              label = "Download example data",
              class = "btn-outline-primary",
              style = "margin-right: 20px;"
              )
          )
        )
      ),
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        upload_plot_countstats_ui(
          id = ns("countStats"),
          title = "Count Stats",
          info.text = "Information about the uploaded counts.",
          caption = "Information about the uploaded counts.",
          height = c("75%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        ),
        upload_plot_phenostats_ui(
          id = ns("phenoStats"),
          title = "Pheno Stats",
          info.text = "Information about the uploaded samples",
          caption = "Information about the uploaded samples.",
          height = c("75%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        ),
        upload_plot_contraststats_ui(
          id = ns("contrastStats"),
          title = "Comparison Stats",
          info.text = "Information about the uploaded comparisons",
          caption = "Information about the uploaded comparisons.",
          height = c("75%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      ),
      bs_alert("In this panel, you can upload your data to the platform. The platform
               requires 3 data files as explained below: a data file containing
               counts/expression (counts.csv), a sample information file (samples.csv)
               and a file specifying the statistical comparisons (comparisons.csv).
               NB Users can now create comparisons from the platform itself, so the
               comparisons.csv file is optional.")
    )   
  )

  comparisons_panel <- shiny::tabPanel(
    "Comparisons",
    bslib::layout_columns(
      col_widths = 12,
      height = "calc(100vh - 200px)",
      heights_equal = "row",
      upload_module_makecontrast_ui(ns("makecontrast")),      
      bs_alert(HTML("Here, you can interactively <b>create comparisons</b> (also called 'contrasts'). Choose a phenotype, then create groups by dragging conditions to the boxes of the 'main' or 'control' group. Give the contrast a name (please keep it short!) and then click 'add comparison'. If you are feeling lucky, you can also try 'auto-comparisons'."))
    )
  )

  batchcorrect_panel <- shiny::tabPanel(
    "BatchEffects",
    bslib::layout_columns(
      col_widths = 12,
      height = "calc(100vh - 180px)",
      heights_equal = "row",
      upload_module_batchcorrect_ui(ns("batchcorrect")),
      bs_alert("Omics data often suffers from batch effect due to experiments done on different days, using different machines or done at different institutes. This will often cause so-called batch effects. Batch correction can clean your data from these 'unwanted variation'. But be careful, batch correction can also be dangerous if not used carefully and can remove valuable real signal. Only adviced for advanced users!")      
    )
  )

  outliers_panel <- shiny::tabPanel(
    "QC/BC",
    bslib::layout_columns(
      col_widths = 12,
      height = "calc(100vh - 180px)",
      heights_equal = "row",
      upload_module_outliers_ui(ns("checkqc")),
      bs_alert("Check for normalization, outliers and batch-effects.")      
    )
  )
  
  compute_panel <- shiny::tabPanel(
    "Compute",
    bs_alert("OK. We now have everything to compute your data. Please name your dataset and give a short description of the experiment. You can select/deselect some computation options but if you do not understand, it is safer to leave the defaults. If you are ready, hit 'Compute'. Computation can take 10-40 minutes depending on the size of your data and number of comparisons."),
    shiny::br(), shiny::br(),
    bslib::layout_columns(
      col_widths = c(1,10,1),
      div(
##        shinyWidgets::prettySwitch(ns("show_batchcorrection"), "Batch correction"),
##        shinyWidgets::prettySwitch(ns("show_checkoutliers"), "Check outliers (beta)")
      ),
      shiny::fillCol(
        height = "100%"
      , style = "padding-right: 40px;",
        upload_module_computepgx_ui(ns("compute"))
      )
    )
  )

  div(
    class = "p-0",
    board_header,
    div(
      style = "position: fixed; right: 0px; width: 160px; margin-top: 10px;",
      shinyWidgets::prettySwitch(ns("expert_mode"), "Expert mode"),
    ),
    shiny::tabsetPanel(
      id = ns("tabs"),
##    upload_select_db,
      upload_panel,
      comparisons_panel,
      outliers_panel,
      batchcorrect_panel,
      compute_panel
    )
  )
}
