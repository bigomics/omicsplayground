##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

AcrossBoard <- function(id, pgx, pgx_dir = reactive(NULL),
                        labeltype = shiny::reactive("feature")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    fullH <- 770
    tabH <- "70vh"

    infotext <- tspan(
      "The <strong>Across Datasets</strong> module enables users to query gene expression
      across multiple datasets stored in a TileDB database. This allows comparing expression
      levels of specific genes across all available experiments without loading individual files.
      <br><br>
      Features:
      <ul>
        <li>Query up to 10 genes simultaneously</li>
        <li>View expression as bar plots colored by dataset</li>
        <li>Compare distributions with boxplots</li>
        <li>Filter samples by phenotype values</li>
        <li>Download data tables with phenotype annotations</li>
      </ul>",
      js = FALSE
    )

    ## ================================================================================
    ## ======================= REACTIVE VALUES ========================================
    ## ================================================================================

    tiledb_path <- reactive({
      req(pgx_dir())
      path <- file.path(pgx_dir(), "counts_tiledb")
      if (dir.exists(path)) return(path)
      parent_path <- file.path(dirname(pgx_dir()), "counts_tiledb")
      if (dir.exists(parent_path)) return(parent_path)
      NULL
    })

    available_datasets <- reactive({
      req(tiledb_path())
      playbase::pgx.listDatasetsTileDB(tiledb_path())
    })

    available_genes <- reactive({
      req(tiledb_path())
      playbase::pgx.listGenesTileDB(tiledb_path())
    })

    phenotypes_by_dataset <- reactive({
      req(tiledb_path())
      playbase::pgx.listPhenotypesByDatasetTileDB(tiledb_path())
    })

    common_phenotypes <- reactive({
      pheno_by_ds <- phenotypes_by_dataset()
      if (length(pheno_by_ds) == 0) return(character(0))

      selected <- input$selected_datasets
      if (length(selected) == 0) selected <- names(pheno_by_ds)

      pheno_lists <- pheno_by_ds[selected]
      if (length(pheno_lists) == 0) return(character(0))

      sort(Reduce(intersect, pheno_lists))
    })

    phenotype_values <- reactive({
      req(input$filter_phenotype)
      req(tiledb_path())
      selected_datasets <- if (length(input$selected_datasets) == 0) NULL else input$selected_datasets
      playbase::pgx.getPhenotypeValuesTileDB(
        tiledb_path(),
        phenotype = input$filter_phenotype,
        datasets = selected_datasets
      )
    })

    query_result <- reactiveVal(NULL)

    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Across Datasets</strong>"),
        shiny::HTML(infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    shiny::observe({
      datasets <- available_datasets()
      if (length(datasets) > 0) {
        shiny::updateSelectizeInput(session, "selected_datasets", choices = datasets, selected = character(0))
      }
    })

    shiny::observe({
      genes <- available_genes()
      if (length(genes) > 0) {
        shiny::updateSelectizeInput(session, "selected_genes", choices = genes, server = TRUE)
      }
    })

    shiny::observe({
      common <- common_phenotypes()
      filter_choices <- if (length(common) > 0) {
        c(setNames("", "(no filter)"), setNames(common, common))
      } else {
        c(setNames("", "(no phenotypes available)"))
      }
      shiny::updateSelectizeInput(session, "filter_phenotype", choices = filter_choices, selected = "")
    })

    shiny::observe({
      filter_pheno <- input$filter_phenotype
      color_choices <- c("dataset" = "dataset")
      if (!is.null(filter_pheno) && filter_pheno != "") {
        color_choices <- c(color_choices, setNames(filter_pheno, filter_pheno))
      }
      shiny::updateSelectizeInput(session, "color_by_phenotype", choices = color_choices, selected = "dataset")
    })

    shiny::observe({
      values <- phenotype_values()
      shiny::updateSelectizeInput(session, "filter_values", choices = values, selected = character(0))
    })

    output$db_info <- shiny::renderText({
      path <- tiledb_path()
      if (is.null(path)) return("No TileDB database found")
      metadata_path <- paste0(path, "_metadata.rds")
      if (!file.exists(metadata_path)) return("Database found (no metadata)")
      metadata <- readRDS(metadata_path)
      n_pheno <- if (!is.null(metadata$n_phenotypes)) metadata$n_phenotypes else 0
      paste0(metadata$n_genes, " genes, ", metadata$n_samples, " samples, ",
             metadata$n_files, " datasets, ", n_pheno, " phenotypes")
    })

    shiny::observeEvent(input$query_button, {
      req(input$selected_genes)
      req(tiledb_path())

      genes <- input$selected_genes
      path <- tiledb_path()
      selected_datasets <- input$selected_datasets
      filter_phenotype <- input$filter_phenotype
      filter_values <- input$filter_values
      color_by <- input$color_by_phenotype

      shiny::withProgress(message = "Querying TileDB...", value = 0.3, {
        counts <- playbase::pgx.queryTileDB(path, genes = genes)
        shiny::incProgress(0.2, message = "Processing data...")

        phenotypes_to_fetch <- character(0)
        if (!is.null(filter_phenotype) && filter_phenotype != "") {
          phenotypes_to_fetch <- c(phenotypes_to_fetch, filter_phenotype)
        }
        if (!is.null(color_by) && color_by != "dataset") {
          phenotypes_to_fetch <- c(phenotypes_to_fetch, color_by)
        }
        phenotypes_to_fetch <- unique(phenotypes_to_fetch)

        df <- if (length(phenotypes_to_fetch) > 0) {
          playbase::pgx.tiledbToPlotDF(counts, tiledb_path = path, phenotypes = phenotypes_to_fetch)
        } else {
          playbase::pgx.tiledbToPlotDF(counts)
        }

        df <- df[!is.na(df$count), ]

        if (length(selected_datasets) > 0) {
          df <- df[df$dataset %in% selected_datasets, ]
        }

        if (!is.null(filter_phenotype) && filter_phenotype != "" &&
            filter_phenotype %in% colnames(df) && length(filter_values) > 0) {
          df <- df[df[[filter_phenotype]] %in% filter_values, ]
        }

        shiny::incProgress(0.3, message = "Done!")
        query_result(df)
      })
    })

    ## ================================================================================
    ## ======================= REACTIVE DATA ==========================================
    ## ================================================================================

    getColorBy <- reactive({
      color_by <- input$color_by_phenotype
      if (is.null(color_by) || color_by == "") "dataset" else color_by
    })

    getPlotData <- reactive({
      df <- query_result()
      if (is.null(df) || nrow(df) == 0) return(NULL)

      if (input$plot_scale == "log2" && any(df$count > 0)) {
        df$count <- log2(df$count + 1)
      }

      color_by <- getColorBy()
      df$color_group <- if (color_by %in% colnames(df)) df[[color_by]] else df$dataset
      attr(df, "color_by") <- color_by
      df
    })

    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    across_plot_barplot_server("barplot", getPlotData = getPlotData, watermark = WATERMARK)
    across_plot_boxplot_server("boxplot", getPlotData = getPlotData, watermark = WATERMARK)
    across_table_data_server("datatable", getPlotData = getPlotData)
  })
}

