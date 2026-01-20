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

    datasets_info_file <- reactive({
      req(pgx_dir())
      path <- file.path(pgx_dir(), "datasets-info.csv")
      if (file.exists(path)) return(path)
      parent_path <- file.path(dirname(pgx_dir()), "datasets-info.csv")
      if (file.exists(parent_path)) return(parent_path)
      NULL
    })

    available_datasets <- reactive({
      req(tiledb_path())
      playbase::pgx.listDatasetsTileDB(tiledb_path())
    })

    dataset_info_table <- reactive({
      req(tiledb_path())
      playbase::pgx.getDatasetInfoTileDB(
        tiledb_path(),
        datasets_info_file = datasets_info_file()
      )
    })

    available_genes <- reactive({
      req(tiledb_path())
      playbase::pgx.listGenesTileDB(tiledb_path())
    })

    phenotypes_by_dataset <- reactive({
      req(tiledb_path())
      playbase::pgx.listPhenotypesByDatasetTileDB(tiledb_path())
    })

    ## Selected datasets (managed via modal)
    selected_datasets <- reactiveVal(character(0))

    common_phenotypes <- reactive({
      pheno_by_ds <- phenotypes_by_dataset()
      if (length(pheno_by_ds) == 0) return(character(0))

      selected <- selected_datasets()
      if (length(selected) == 0) selected <- names(pheno_by_ds)

      pheno_lists <- pheno_by_ds[selected]
      if (length(pheno_lists) == 0) return(character(0))

      sort(Reduce(intersect, pheno_lists))
    })

    phenotype_values <- reactive({
      req(input$filter_phenotype)
      req(tiledb_path())
      sel_datasets <- if (length(selected_datasets()) == 0) NULL else selected_datasets()
      playbase::pgx.getPhenotypeValuesTileDB(
        tiledb_path(),
        phenotype = input$filter_phenotype,
        datasets = sel_datasets
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

    ## ================================================================================
    ## ======================= DATASET SELECTOR MODAL =================================
    ## ================================================================================

    ## Button label
    output$dataset_button_label <- shiny::renderUI({
      sel <- selected_datasets()
      total <- length(available_datasets())
      if (length(sel) == 0) {
        shiny::span("All datasets", style = "color: #666;")
      } else {
        shiny::span(paste0(length(sel), " of ", total, " selected"))
      }
    })

    ## Selected datasets text (truncated list of names)
    output$selected_datasets_text <- shiny::renderUI({
      sel <- selected_datasets()
      if (length(sel) == 0) {
        return(shiny::tags$em("(querying all datasets)"))
      }

      ## Build truncated string with max ~60 characters
      max_chars <- 60
      result <- ""
      shown_count <- 0

      for (name in sel) {
        if (shown_count == 0) {
          candidate <- name
        } else {
          candidate <- paste0(result, ", ", name)
        }

        if (nchar(candidate) > max_chars && shown_count > 0) {
          break
        }

        result <- candidate
        shown_count <- shown_count + 1
      }

      remaining <- length(sel) - shown_count
      if (remaining > 0) {
        result <- paste0(result, " +", remaining, " more")
      }

      shiny::tags$span(result)
    })

    ## Store selection state before opening modal (for Cancel to revert)
    selection_before_modal <- reactiveVal(character(0))

    ## Open modal (show persistent div)
    shiny::observeEvent(input$open_dataset_modal, {
      ## Store current selection for Cancel
      selection_before_modal(selected_datasets())
      ## Show the modal
      shinyjs::show("dataset_modal_container")
    })

    ## Helper to hide modal
    hideModal <- function() {
      shinyjs::hide("dataset_modal_container")
    }

    ## Modal selection count - reads directly from DT selection (no re-render)
    output$modal_selection_count <- shiny::renderText({
      sel_rows <- input$dataset_table_rows_selected
      total <- length(available_datasets())
      n_selected <- length(sel_rows)
      if (n_selected == 0) {
        paste0("No datasets selected (will query all ", total, " datasets)")
      } else {
        paste0(n_selected, " of ", total, " datasets selected")
      }
    })

    ## Dataset table in modal - persistent, only re-renders when data changes
    output$dataset_table <- DT::renderDataTable({
      df <- dataset_info_table()
      if (is.null(df) || nrow(df) == 0) {
        return(DT::datatable(data.frame(Message = "No datasets found")))
      }

      ## Select columns to display
      display_cols <- c("dataset", "nsamples", "datatype", "organism", "description")
      display_cols <- intersect(display_cols, colnames(df))
      df <- df[, display_cols, drop = FALSE]

      ## Truncate long descriptions
      if ("description" %in% colnames(df)) {
        df$description <- sapply(df$description, function(x) {
          if (is.na(x)) return("")
          if (nchar(x) > 80) paste0(substr(x, 1, 77), "...") else x
        })
      }

      ## Convert categorical columns to factors for dropdown filters
      if ("datatype" %in% colnames(df)) {
        df$datatype <- factor(df$datatype)
      }
      if ("organism" %in% colnames(df)) {
        df$organism <- factor(df$organism)
      }

      ## Initial selection from current state (only affects first render)
      init_sel <- isolate(selected_datasets())
      sel_rows <- which(df$dataset %in% init_sel)

      DT::datatable(
        df,
        class = "compact hover",
        rownames = FALSE,
        width = "100%",
        selection = list(mode = "multiple", selected = sel_rows),
        filter = list(position = "top", clear = TRUE),
        options = list(
          dom = "tip",
          pageLength = 15,
          scrollX = FALSE,
          autoWidth = TRUE
        )
      )
    })

    ## Select all visible
    shiny::observeEvent(input$modal_select_all, {
      df <- dataset_info_table()
      if (!is.null(df) && nrow(df) > 0) {
        ## Get filtered rows if filter is active
        visible_rows <- input$dataset_table_rows_all
        if (is.null(visible_rows)) visible_rows <- seq_len(nrow(df))
        DT::dataTableProxy("dataset_table") |> DT::selectRows(visible_rows)
      }
    })

    ## Clear selection
    shiny::observeEvent(input$modal_clear, {
      DT::dataTableProxy("dataset_table") |> DT::selectRows(NULL)
    })

    ## Apply selection - save current DT selection and close
    shiny::observeEvent(input$modal_apply, {
      df <- dataset_info_table()
      sel_rows <- input$dataset_table_rows_selected
      if (!is.null(df) && length(sel_rows) > 0) {
        selected_datasets(df$dataset[sel_rows])
      } else {
        selected_datasets(character(0))
      }
      hideModal()
    })

    ## Cancel - revert to previous selection and close
    shiny::observeEvent(input$modal_cancel, {
      ## Restore previous selection in DT
      df <- dataset_info_table()
      prev_sel <- selection_before_modal()
      if (!is.null(df) && nrow(df) > 0) {
        sel_rows <- which(df$dataset %in% prev_sel)
        DT::dataTableProxy("dataset_table") |> DT::selectRows(sel_rows)
      }
      hideModal()
    })

    ## Cancel button (second one in footer)
    shiny::observeEvent(input$modal_cancel2, {
      df <- dataset_info_table()
      prev_sel <- selection_before_modal()
      if (!is.null(df) && nrow(df) > 0) {
        sel_rows <- which(df$dataset %in% prev_sel)
        DT::dataTableProxy("dataset_table") |> DT::selectRows(sel_rows)
      }
      hideModal()
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
      sel_datasets <- selected_datasets()
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

        if (length(sel_datasets) > 0) {
          df <- df[df$dataset %in% sel_datasets, ]
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

