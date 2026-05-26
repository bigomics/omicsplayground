##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

downloadButton2 <- function(outputId, label = "Download", class = NULL, ...) {
  aTag <- shiny::tags$a(
    id = outputId,
    class = paste("btn btn-default shiny-download-link", class),
    href = "", target = "_blank", download = NA,
    shiny::icon("file-csv"), label, ...
  )
}

LoadingInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings()
}

LoadingUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  board_header <- fillRow(
    flex = c(NA, 1, NA),
    shiny::div(
      id = "navheader-current-section",
      HTML("Data Library &nbsp;"),
      shiny::actionLink(
        ns("module_info"), "",
        icon = shiny::icon("youtube"),
        style = "color: #ccc;"
      )
    ),
    div(),
    shiny::div(shiny::uiOutput(ns("pgx_stats_ui")), id = "navheader-dataset-stats")
  )

  user_tabpanel <- shiny::tabPanel(
    "My Datasets",
    bslib::layout_columns(
      col_widths = 12,
      height = "calc(100vh - 141px)",
      uiOutput(ns("sharing_alert")),
      div(
        shiny::actionButton(
          ns("newuploadbutton"),
          label = "Upload new",
          icon = icon("upload"),          
          class = "btn btn-outline-primary"
        ),        
        shiny::actionButton(
          ns("loadbutton"),
          label = "Load selected",
          icon = icon("file-import"),
          class = "btn btn-primary"
        )
      ),
      bslib::layout_columns(
        col_widths = c(8, 4),
        loading_table_datasets_ui(
          ns("pgxtable"),
          title = "Available datasets",
          info.text = "This table contains information about all available datasets within the platform. For each dataset, it reports a brief description as well as the total number of samples, genes, genesets, corresponding phenotypes and the creation date.",
          caption = "Table of datasets available in the platform.",
          height = c("100vh", 700),
          width = c("100%", "100%")
        ),
        loading_tsne_ui(
          ns("tsne"),
          title = "Signature t-SNE",
          info.text = "Scatter plot displaying the t-SNE clustering of the available contrasts. Each dot corresponds to a specific comparison.",
          info.methods = "t-SNE is a non-linear dimensionality reduction method that enables visualization of high-dimensional data in a low-dimensional space, typically 2D or 3D. Unlike linear dimensionality reduction techniques like PCA, t-SNE may separate data that is not linearly separable. Signatures/datasets that are clustered closer together, are more similar. Performed using the Rtsne R package [1].",
          info.references = list(
            list(
              "Krijthe JH (2015). Rtsne: T-Distributed Stochastic Neighbor Embedding using Barnes-Hut Implementation. R package version 0.17",
              "https://doi.org/10.32614/CRAN.package.Rtsne"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
          caption = "Similarity clustering of fold-change signatures colored by data sets using t-SNE.",
          height = c("100vh", "70vh"),
          width = c("auto", "100%")
        )
      ) ## end of 7fr-5fr
    )
  )

  public_tabpanel <- shiny::tabPanel(
    opt$PUBLIC_DATASETS_LABEL,
    bslib::layout_columns(
      col_widths = 12,
      height = "calc(100vh - 181px)",
      bs_alert(
        if (opt$ENABLE_PUBLIC_LOAD) {
          paste0("This panel shows all <b>", tolower(opt$PUBLIC_DATASETS_LABEL), "</b>. You can select a public dataset and click <b>Load selected</b> to load it directly for analysis (without importing), or click <b>Import Dataset</b> to copy it to your library. The <b>Signature t-SNE</b> shows similarity clustering of fold-change signatures using t-SNE.")
        } else {
          paste0("This panel shows all <b>", tolower(opt$PUBLIC_DATASETS_LABEL), "</b>. You can select a public dataset and click <b>Import Dataset</b> to copy that dataset to your library for further analysis. The <b>Signature t-SNE</b> shows similarity clustering of fold-change signatures using t-SNE.")
        },
        translate = FALSE, html = TRUE
      ),
      bslib::layout_columns(
        col_widths = c(8, 4),
        height = "calc(100vh - 181px)",
        loading_table_datasets_public_ui(
          ns("pgxtable_public"),
          title = opt$PUBLIC_DATASETS_LABEL,
          info.text = paste0("This table shows available ", tolower(opt$PUBLIC_DATASETS_LABEL), " within the platform. For each dataset, it reports a brief description as well as the total number of samples, genes, gene sets (or pathways), corresponding phenotypes and the creation date."),
          caption = paste0("Table with ", tolower(opt$PUBLIC_DATASETS_LABEL), " available in the platform."),
          ## height = c("calc(100vh - 330px)", 700),
          height = c("100%", 700),
          width = c("100%", "100%"),
          load_button = opt$ENABLE_PUBLIC_LOAD,
          delete_button = TRUE
        ),
        loading_tsne_ui(
          ns("tsne_public"),
          title = "Signature t-SNE",
          info.text = "Scatter plot displaying the t-SNE clustering of the available contrasts. Each dot corresponds to a specific comparison.",
          info.methods = "t-SNE is a non-linear dimensionality reduction method that enables visualization of high-dimensional data in a low-dimensional space, typically 2D or 3D. Unlike linear dimensionality reduction techniques like PCA, t-SNE may separate data that is not linearly separable. Signatures/datasets that are clustered closer together, are more similar. Performed using the Rtsne R package [1].",
          info.references = list(
            list(
              "Krijthe JH (2015). Rtsne: T-Distributed Stochastic Neighbor Embedding using Barnes-Hut Implementation. R package version 0.17",
              "https://doi.org/10.32614/CRAN.package.Rtsne"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
          caption = "Similarity clustering of fold-change signatures colored by data sets using t-SNE.",
          ## height = c("calc(100vh - 330px)", 700),
          height = c("100%", 700),
          width = c("auto", "100%")
        )
      ) ## end of 7fr-5fr
    )
  ) ## end of Public tabPanel

  archive_tabpanel <- shiny::tabPanel(
    "Data archive",
    value = "archive_tab",
    bslib::layout_columns(
      col_widths = 12,
      height = "calc(100vh - 181px)",
      bs_alert(
        if (opt$ENABLE_PUBLIC_LOAD) {
          "This panel shows all <b>Archived datasets</b>. You can select an archived dataset and click <b>Load selected</b> to load it directly for analysis (without importing), or click <b>Import Dataset</b> to copy it to your library. The <b>Signature t-SNE</b> shows similarity clustering of fold-change signatures using t-SNE."
        } else {
          "This panel shows all <b>Archived datasets</b>. You can select an archived dataset and click <b>Import Dataset</b> to copy that dataset to your library for further analysis. The <b>Signature t-SNE</b> shows similarity clustering of fold-change signatures using t-SNE."
        },
        translate = FALSE, html = TRUE
      ),
      bslib::layout_columns(
        col_widths = c(8, 4),
        height = "calc(100vh - 181px)",
        loading_table_datasets_public_ui(
          ns("pgxtable_archive"),
          title = "Archived datasets",
          info.text = "This table shows available datasets within the platform that have been archived. For each dataset, it reports a brief description as well as the total number of samples, genes, gene sets (or pathways), corresponding phenotypes and the creation date.",
          caption = "Table with archived datasets available in the platform.",
          ## height = c("calc(100vh - 330px)", 700),
          height = c("100%", 700),
          width = c("100%", "100%"),
          delete_button = TRUE,
          load_button = opt$ENABLE_PUBLIC_LOAD
        ),
        loading_tsne_ui(
          ns("tsne_archive"),
          title = "Signature t-SNE",
          info.text = "Scatter plot displaying the t-SNE clustering of the available contrasts. Each dot corresponds to a specific comparison.",
          info.methods = "t-SNE is a non-linear dimensionality reduction method that enables visualization of high-dimensional data in a low-dimensional space, typically 2D or 3D. Unlike linear dimensionality reduction techniques like PCA, t-SNE may separate data that is not linearly separable. Signatures/datasets that are clustered closer together, are more similar. Performed using the Rtsne R package [1].",
          info.references = list(
            list(
              "Krijthe JH (2015). Rtsne: T-Distributed Stochastic Neighbor Embedding using Barnes-Hut Implementation. R package version 0.17",
              "https://doi.org/10.32614/CRAN.package.Rtsne"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
          caption = "Similarity clustering of fold-change signatures colored by data sets using t-SNE.",
          ## height = c("calc(100vh - 330px)", 700),
          height = c("100%", 700),
          width = c("auto", "100%")
        )
      ) ## end of 7fr-5fr
    ),
    info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/data/#archived-datasets",
    caption = "Table with archived datasets available in the platform.",
    height = "calc(100vh - 181px)",
    width = "100%"
  )

  sharing_tabpanel <- shiny::tabPanel(
    "Shared datasets",
    value = "sharing_tab",
    bslib::layout_columns(
      col_widths = 12,
      #height = "calc(100vh - 181px)",
      height = "100%",
      row_heights = c("auto",1),
      bs_alert("This Sharing panel shows <strong>received datasets</strong> that are not yet imported to your library, and your <strong>shared datasets</strong> that are still waiting to be accepted by the receiver. Please accept or refuse each received file, and/or resend a message or cancel your shared datasets using the action buttons on the right of the tables."),
      bslib::layout_columns(
        col_widths = c(6,6),
        height = "100%",
        upload_module_received_ui(ns("received")),
        upload_module_shared_ui(ns("shared"))        
      )
    )
  )
  
  ## ------------------------------------------------------------------------

  ## disable/hide public tabpanel if public folder does not exists
  public_dir <- file.path(OPG, "data_public")
  if (!dir.exists(public_dir)) {
    public_tabpanel <- NULL
  }

  ## ============================ Board object ===========================
  div(
    class = "p-0",
    board_header,
    shiny::tabsetPanel(
      id = ns("tabs"),
      user_tabpanel,
      public_tabpanel,
      archive_tabpanel,
      sharing_tabpanel      
    )
  )
}

## ====================================================================
## ====================================================================
## ====================================================================

