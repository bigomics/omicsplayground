##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

AcrossInputs <- function(id) {

  ns <- shiny::NS(id) ## namespace
  ## NOTE: rendered inside the top-level "Across datasets" panel (no Dashboard
  ## settings sidebar), so we use a plain container instead of bigdash::tabSettings().
  shiny::div(
    class = "across-inputs",
    ## Dataset selector button (opens modal)
    shiny::tags$label("Select datasets:", class = "control-label"),
    shiny::div(
      style = "margin-bottom: 5px;",
      shiny::actionButton(
        ns("open_dataset_modal"),
        shiny::uiOutput(ns("dataset_button_label"), inline = TRUE),
        class = "btn-outline-secondary btn-sm",
        style = "width: 100%; text-align: left;",
        icon = shiny::icon("database")
      )
    ),
    shiny::div(
      class = "across-selected-datasets-text",
      shiny::uiOutput(ns("selected_datasets_text"))
    ),
    shiny::br(),
    withTooltip(
      shiny::selectizeInput(
        ns("selected_genes"),
        "Select genes:",
        choices = NULL,
        multiple = TRUE,
        options = list(
          maxItems = 10,
          placeholder = "Type gene names..."
        )
      ),
      "Select up to 10 genes to query across selected datasets.",
      placement = "right", options = list(container = "body")
    ),
    shiny::br(),
    withTooltip(
      shiny::actionButton(
        ns("query_button"),
        label = "Query",
        class = "btn-outline-primary",
        icon = icon("search")
      ),
      "Click to query selected genes across selected datasets.",
      placement = "right"
    ),
    shiny::br(),
    shiny::br(),
    bslib::accordion(
      id = ns("across_accordion"),
      open = FALSE,
      bslib::accordion_panel(
        "Split by Phenotype",
        icon = icon("th-list", lib = "glyphicon"),
        shiny::tags$small(
          shiny::tags$em("Expression is always shown per dataset. Optionally split each dataset into sub-groups by a phenotype.")
        ),
        shiny::br(),
        shiny::br(),
        withTooltip(
          shiny::radioButtons(
            ns("pheno_mode"),
            "Phenotypes to offer:",
            choices = c("All datasets" = "union", "Shared only" = "intersection"),
            selected = "union",
            inline = TRUE
          ),
          paste(
            "All: every phenotype found in any selected dataset (datasets missing it form an 'n/a' group);",
            "several can be combined. Shared: only phenotypes present in all selected datasets."
          ),
          placement = "right", options = list(container = "body")
        ),
        ## Selector is rendered server-side: single-select in 'Shared' mode,
        ## multi-select in 'All' mode (see output$split_by_ui).
        shiny::uiOutput(ns("split_by_ui"))
      ),
      bslib::accordion_panel(
        "Options",
        icon = icon("cog", lib = "glyphicon"),
        shiny::div(
          style = "display: flex; align-items: center; gap: 6px;",
          withTooltip(
            shiny::radioButtons(
              ns("value_type"),
              "Value type:",
              choices = c("Counts" = "count", "Z-score" = "zscore"),
              selected = "count",
              inline = TRUE
            ),
            "Select value type: raw counts or z-scores (per-dataset standardized values).",
            placement = "right", options = list(container = "body")
          ),
          withTooltip(
            shiny::actionLink(
              ns("zscore_info"),
              label = NULL,
              icon = shiny::icon("circle-info"),
              style = "color: #004ca7;"
            ),
            "What is a z-score, and why use it across datasets?",
            placement = "top"
          )
        ),
        withTooltip(
          shiny::radioButtons(
            ns("plot_scale"),
            "Scale:",
            choices = c("linear", "log2"),
            selected = "linear",
            inline = TRUE
          ),
          "Select scale for expression values.",
          placement = "right", options = list(container = "body")
        )
      )
    ),
    shiny::br(),
    shiny::hr(),
    shiny::tags$small(
      shiny::tags$strong("Database info:"),
      shiny::br(),
      shiny::textOutput(ns("db_info"), inline = TRUE)
    )
  )
}

AcrossUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  ## Tab content fills the column height (see the .across-content rules in
  ## scss/components/_across.scss) so the plot card's bottom lines up with the
  ## inputs card on the left, independent of header/tab-strip pixel heights.
  fullH <- "100%"
  tabH <- "70vh"

  tabs <- shiny::tabsetPanel(
    id = ns("tabs1"),
    shiny::tabPanel(
      "Expression by Sample",
      bslib::layout_columns(
        col_widths = 12,
        height = fullH,
        bs_alert("View gene expression values across all samples from multiple datasets stored in the TileDB database."),
        bslib::layout_columns(
          col_widths = 12,
          height = fullH,
          across_plot_barplot_ui(
            id = ns("barplot"),
            title = "Expression by Sample",
            info.text = "Bar plot showing expression values for selected genes across all samples. Samples are colored by their source dataset.",
            info.methods = "Expression values are retrieved from a TileDB database containing counts from multiple PGX files. Values are displayed as raw counts or log2-transformed based on the scale setting.",
            info.extra_link = NULL,
            width = c("auto", "100%"),
            height = c("100%", "75vh")
          )
        )
      )
    ),
    shiny::tabPanel(
      "Expression by Dataset",
      bslib::layout_columns(
        col_widths = 12,
        height = fullH,
        bs_alert("Compare gene expression distributions across different datasets using boxplots."),
        bslib::layout_columns(
          col_widths = 12,
          height = fullH,
          across_plot_boxplot_ui(
            id = ns("boxplot"),
            title = "Distribution by Dataset",
            info.text = "Boxplot showing the distribution of expression values for selected genes across different datasets.",
            info.methods = "Expression values are grouped by source dataset and displayed as boxplots. Jittered points can be overlaid to show individual sample values.",
            info.extra_link = NULL,
            width = c("auto", "100%"),
            height = c("100%", "75vh")
          )
        )
      )
    ),
    shiny::tabPanel(
      "Data Table",
      bslib::layout_columns(
        col_widths = 12,
        height = fullH,
        bs_alert("View and download the raw expression data for selected genes."),
        bslib::layout_columns(
          col_widths = 12,
          height = fullH,
        across_table_data_ui(
          id = ns("datatable"),
          title = "Expression Data",
            height = c("100%", "75vh"),
            width = c("auto", "100%")
          )
        )
      )
    )
  )

  content <- div(
    class = "across-content",
    tabs,
    ## Persistent dataset selector modal (hidden by default)
    shiny::tags$div(
      id = ns("dataset_modal_container"),
      style = "display: none;",
      ## NOTE: do NOT use the Bootstrap "fade show" classes here. Bootstrap's
      ## modal data-api keys off `.modal.show`: on every click of a
      ## `data-bs-toggle="modal"` trigger (e.g. the plot fullscreen/zoom
      ## buttons) it does `Modal.getInstance(findOne('.modal.show')).hide()`.
      ## This hand-rolled overlay has no Bootstrap instance, so a persistent
      ## `.modal.show` made that call throw and silently broke every fullscreen
      ## button. `modal` alone gives the same static "shown" look (the overlay's
      ## visibility is controlled by toggling the container's display).
      shiny::tags$div(
        class = "modal across-dataset-modal-overlay",
        style = "display: block;",
        tabindex = "-1",
        shiny::tags$div(
          class = "modal-dialog modal-xl across-dataset-modal-dialog",
          shiny::tags$div(
            class = "modal-content",
            shiny::tags$div(
              class = "modal-header",
              shiny::tags$h5(class = "modal-title", "Select Datasets"),
              shiny::actionButton(
                ns("modal_cancel"),
                label = NULL,
                icon = shiny::icon("times"),
                class = "btn-close",
                style = "background: none; border: none;"
              )
            ),
            shiny::tags$div(
              class = "modal-body",
              shiny::div(
                style = "margin-bottom: 10px;",
                shiny::textOutput(ns("modal_selection_count"))
              ),
              DT::dataTableOutput(ns("dataset_table"))
            ),
            shiny::tags$div(
              class = "modal-footer",
              shiny::div(
                style = "margin-right: auto;",
                shiny::actionButton(ns("modal_select_all"), "Select All", class = "btn btn-sm btn-outline-secondary"),
                shiny::actionButton(ns("modal_clear"), "Clear", class = "btn btn-sm btn-outline-secondary")
              ),
              shiny::actionButton(ns("modal_cancel2"), "Cancel", class = "btn btn-secondary"),
              shiny::actionButton(ns("modal_apply"), "Apply", class = "btn btn-primary")
            )
          )
        )
      )
    )
  )

  ## Board header spans the full width; the inputs (left) and content (right)
  ## align on the same baseline below it.
  bslib::layout_columns(
    col_widths = 12,
    gap = 0,
    height = "calc(100vh - 60px)",
    row_heights = list("auto", 1),
    boardHeader(title = "Across Datasets", info_link = ns("info")),
    bslib::layout_columns(
      col_widths = c(3, 9),
      fill = TRUE,
      bslib::card(
        bslib::card_body(
          AcrossInputs(id)
        )
      ),
      content
    )
  )
}

