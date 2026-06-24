##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

AcrossBoard <- function(id, pgx, pgx_dir = reactive(NULL),
                        labeltype = shiny::reactive("feature"),
                        tiledb_refresh = shiny::reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    fullH <- 770
    tabH <- "70vh"

    infotext <- tspan(
      "The <strong>Across Datasets</strong> module lets you compare feature values
      across all of your datasets at once, without loading each one individually. Compare
      the levels of specific features across every available experiment side by side.
      <br><br>
      What you can do:
      <ul>
        <li>Query up to 10 features at once</li>
        <li>View values as bar plots colored by dataset</li>
        <li>Compare distributions with boxplots</li>
        <li>Split each dataset into sub-groups by phenotype</li>
        <li>Download data tables with phenotype annotations</li>
      </ul>",
      js = FALSE
    )

    ## ================================================================================
    ## ======================= REACTIVE VALUES ========================================
    ## ================================================================================

    tiledb_path <- reactive({
      ## Re-check the filesystem when the background TileDB build completes:
      ## without this dependency the first (empty) result is cached and the
      ## board stays blank until the user logs out and back in.
      tiledb_refresh()
      req(pgx_dir())
      path <- file.path(pgx_dir(), "counts_tiledb")
      if (dir.exists(path)) return(path)
      parent_path <- file.path(dirname(pgx_dir()), "counts_tiledb")
      if (dir.exists(parent_path)) return(parent_path)
      NULL
    })

    datasets_info_file <- reactive({
      tiledb_refresh()
      req(pgx_dir())
      path <- file.path(pgx_dir(), "datasets-info.csv")
      if (file.exists(path)) return(path)
      parent_path <- file.path(dirname(pgx_dir()), "datasets-info.csv")
      if (file.exists(parent_path)) return(parent_path)
      NULL
    })

    available_datasets <- reactive({
      req(tiledb_path())
      playbase::listDatasetsTileDB(tiledb_path())
    })

    dataset_info_table <- reactive({
      req(tiledb_path())
      playbase::getDatasetInfoTileDB(
        tiledb_path(),
        datasets_info_file = datasets_info_file()
      )
    })

    available_genes <- reactive({
      req(tiledb_path())
      playbase::listGenesTileDB(tiledb_path())
    })

    phenotypes_by_dataset <- reactive({
      req(tiledb_path())
      playbase::listPhenotypesByDatasetTileDB(tiledb_path())
    })

    ## Selected datasets (managed via modal)
    selected_datasets <- reactiveVal(character(0))

    ## Phenotypes offered for splitting. In "intersection" mode only phenotypes
    ## present in every selected dataset are offered; in "union" mode every
    ## phenotype found in any selected dataset is offered (datasets missing it
    ## form an "n/a" group at plot time).
    candidate_phenotypes <- reactive({
      pheno_by_ds <- phenotypes_by_dataset()
      if (length(pheno_by_ds) == 0) return(character(0))

      selected <- selected_datasets()
      if (length(selected) == 0) selected <- names(pheno_by_ds)

      pheno_lists <- pheno_by_ds[selected]
      pheno_lists <- pheno_lists[!vapply(pheno_lists, is.null, logical(1))]
      if (length(pheno_lists) == 0) return(character(0))

      mode <- if (is.null(input$pheno_mode)) "union" else input$pheno_mode
      if (mode == "intersection") {
        sort(Reduce(intersect, pheno_lists))
      } else {
        sort(unique(unlist(pheno_lists)))
      }
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

    shiny::observeEvent(input$zscore_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Z-score: comparing features across datasets</strong>"),
        shiny::HTML(
          "<p><b>What is a z-score?</b><br>
          A z-score standardizes each feature's values <i>within a dataset</i>: for every
          feature we subtract that feature's mean across the dataset's samples and divide by its
          standard deviation. The result has mean&nbsp;0 and standard deviation&nbsp;1, so a
          z-score of <b>+2</b> means &ldquo;two standard deviations above this feature's average
          in that dataset&rdquo; and <b>&minus;1</b> means &ldquo;one SD below average&rdquo;,
          independent of the original units.</p>
          <p><b>Why use it to compare across datasets?</b><br>
          Raw counts are not directly comparable between datasets &mdash; different
          experiments differ in sequencing depth, platform, normalization and dynamic range,
          so the same absolute value can mean very different things. Standardizing per dataset
          removes these dataset-specific scale and batch effects and puts every dataset on the
          same unit-free scale. That lets you compare a feature's <i>relative</i> expression
          pattern (high vs. low, and how variable) on equal footing across experiments &mdash;
          which is exactly what across-datasets analysis needs.</p>
          <p>Use <b>Counts</b> when you care about absolute levels within a single dataset;
          use <b>Z-score</b> when comparing patterns across datasets.</p>"
        ),
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

      ## Select columns to display (base + metadata columns)
      base_cols <- c("dataset", "nsamples", "datatype", "organism", "description")
      metadata_cols <- grep("^metadata_", colnames(df), value = TRUE)
      display_cols <- c(base_cols, metadata_cols)
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
      ## Ensure proper character conversion and NA handling
      categorical_cols <- c("datatype", "organism", metadata_cols)
      for (col in categorical_cols) {
        if (col %in% colnames(df)) {
          vals <- as.character(df[[col]])
          vals[is.na(vals) | vals == "" | vals == "NA"] <- NA
          unique_vals <- unique(vals[!is.na(vals)])
          if (length(unique_vals) > 0) {
            df[[col]] <- factor(vals, levels = unique_vals)
          } else {
            df[[col]] <- factor(vals)
          }
        }
      }

      ## Clean up metadata column names for display (remove "metadata_" prefix)
      clean_names <- colnames(df)
      clean_names <- gsub("^metadata_", "", clean_names)
      clean_names <- gsub("_", " ", clean_names)
      colnames(df) <- clean_names

      ## Initial selection from current state (only affects first render)
      init_sel <- isolate(selected_datasets())
      sel_rows <- which(df$dataset %in% init_sel)

      ## Prepend an empty column that is rendered as a selection checkbox (see
      ## .across-select-col in _across.scss). Selection is still driven by DT's
      ## native row selection, so clicking the checkbox toggles the row and all
      ## the existing selection logic (count, Select All / Clear, Apply) keeps
      ## working unchanged.
      ## Note: the column name must be non-empty. DT's server-side filter maps
      ## request columns to data columns by name, and an empty name collapses to
      ## index 0, triggering "argument is of length 0". A single space renders as
      ## a blank header while still being a valid, matchable name.
      df <- cbind(data.frame(.sel = rep("", nrow(df)), stringsAsFactors = FALSE), df)
      colnames(df)[1] <- " "

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
          scrollX = TRUE,
          autoWidth = TRUE,
          columnDefs = list(
            list(
              targets = 0,
              orderable = FALSE,
              searchable = FALSE,
              className = "across-select-col"
            )
          )
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

    ## Split-by-phenotype selector. Single-select in "Shared" (intersection)
    ## mode; multi-select in "All" (union) mode, where several phenotypes can be
    ## combined per dataset. Leaving it empty means one group per dataset.
    output$split_by_ui <- shiny::renderUI({
      mode <- if (is.null(input$pheno_mode)) "union" else input$pheno_mode
      phenos <- candidate_phenotypes()
      if (mode == "intersection") {
        ## Shared mode: a true single-select with an explicit "Ungrouped" entry
        ## as the default, so the user can switch back to no grouping just by
        ## picking it (no backspacing). The "__none__" sentinel is not a real
        ## column name, so getPlotData treats it as "no split".
        input_el <- shiny::selectizeInput(
          ns("split_by"),
          "Split by phenotype:",
          choices = c("Ungrouped (one group per dataset)" = "__none__", setNames(phenos, phenos)),
          selected = "__none__",
          multiple = FALSE
        )
      } else {
        ## All mode: multi-select; an empty selection means ungrouped. Pick one
        ## or more phenotypes to combine.
        input_el <- shiny::selectizeInput(
          ns("split_by"),
          "Split by phenotype:",
          choices = setNames(phenos, phenos),
          selected = character(0),
          multiple = TRUE,
          options = list(placeholder = "Ungrouped (one group per dataset) - pick one or more")
        )
      }
      withTooltip(
        input_el,
        "Split each dataset into sub-groups by the selected phenotype(s). Choose 'Ungrouped' for one group per dataset.",
        placement = "right", options = list(container = "body")
      )
    })

    output$db_info <- shiny::renderText({
      path <- tiledb_path()
      if (is.null(path)) return("No datasets found")
      metadata_path <- paste0(path, "_metadata.rds")
      if (!file.exists(metadata_path)) return("Database found (no metadata)")
      metadata <- readRDS(metadata_path)
      n_pheno <- if (!is.null(metadata$n_phenotypes)) metadata$n_phenotypes else 0
      paste0(metadata$n_genes, " features, ", metadata$n_samples, " samples, ",
             metadata$n_files, " datasets, ", n_pheno, " phenotypes")
    })

    ## Helper function to run the TileDB query
    runTileDBQuery <- function() {
      genes <- input$selected_genes
      path <- tiledb_path()
      sel_datasets <- selected_datasets()
      value_type <- if (is.null(input$value_type)) "count" else input$value_type

      ## Fetch every phenotype available in the (selected) datasets up front, so
      ## that changing the "Split by" selection only recomputes groups in
      ## getPlotData and does NOT require another TileDB query. The expensive
      ## part is the counts query; phenotypes come from in-memory metadata.
      pheno_by_ds <- phenotypes_by_dataset()
      ds_for_pheno <- if (length(sel_datasets) > 0) sel_datasets else names(pheno_by_ds)
      all_phenotypes <- sort(unique(unlist(pheno_by_ds[ds_for_pheno])))

      shiny::withProgress(message = "Querying your datasets...", value = 0.3, {
        ## Push the dataset selection down to TileDB: query only the samples that
        ## belong to the selected datasets (samples are "<dataset>::<sample>").
        ## Falls back to all samples (gene push-down only) if unavailable.
        sample_filter <- NULL
        if (length(sel_datasets) > 0) {
          all_s <- tryCatch(playbase::listSamplesTileDB(path), error = function(e) NULL)
          if (!is.null(all_s)) {
            sample_filter <- all_s[playbase::getDatasetFromSample(all_s) %in% sel_datasets]
          }
        }
        counts <- playbase::queryTileDB(path, genes = genes, samples = sample_filter, value = value_type)
        shiny::incProgress(0.2, message = "Processing data...")

        df <- if (length(all_phenotypes) > 0) {
          playbase::tiledbToPlotDF(counts, tiledb_path = path, phenotypes = all_phenotypes)
        } else {
          playbase::tiledbToPlotDF(counts)
        }

        df <- df[!is.na(df$count), ]

        if (length(sel_datasets) > 0) {
          df <- df[df$dataset %in% sel_datasets, ]
        }

        shiny::incProgress(0.3, message = "Done!")
        attr(df, "value_type") <- value_type
        query_result(df)
      })
    }

    shiny::observeEvent(input$query_button, {
      req(input$selected_genes)
      req(tiledb_path())
      runTileDBQuery()
    })

    ## Re-query when value_type changes (only if a query has already been made)
    shiny::observeEvent(input$value_type, {
      req(query_result())  # Only re-query if we already have results
      req(input$selected_genes)
      req(tiledb_path())
      runTileDBQuery()
    }, ignoreInit = TRUE)

    ## Disable scale selector when z-score is selected (z-scores should not be log-transformed)
    shiny::observe({
      if (isTRUE(input$value_type == "zscore")) {
        shinyjs::disable("plot_scale")
      } else {
        shinyjs::enable("plot_scale")
      }
    })

    ## ================================================================================
    ## ======================= REACTIVE DATA ==========================================
    ## ================================================================================

    getPlotData <- reactive({
      df <- query_result()
      if (is.null(df) || nrow(df) == 0) return(NULL)

      value_type <- if (is.null(input$value_type)) "count" else input$value_type
      plot_scale <- ifelse(value_type == "zscore", "linear", input$plot_scale)

      if (plot_scale == "log2" && any(df$count > 0)) {
        df$count <- log2(df$count + 1)
      }

      ## Compute the split here (not in the query) so changing "Split by" or the
      ## Shared/All toggle updates the plots without re-querying TileDB. All
      ## phenotypes were fetched up front; keep only the selected one(s) so the
      ## data table is not cluttered with every phenotype column.
      split_by <- input$split_by
      split_by <- split_by[!is.na(split_by) & nzchar(split_by)]
      present <- intersect(split_by, colnames(df))
      has_split <- length(present) > 0

      std_cols <- c("gene", "sample", "sample_short", "count", "dataset")
      df <- df[, c(intersect(std_cols, colnames(df)), present), drop = FALSE]

      ## Plots always group by dataset. When a split is active the sub-group
      ## (and colour) is the phenotype split; otherwise it is the dataset
      ## itself, giving one box/bar group per dataset. For several phenotypes
      ## (union mode) the group is the combination of the non-missing values a
      ## sample has; samples lacking all of them fall into an "n/a" group.
      if (has_split) {
        parts <- lapply(present, function(p) {
          v <- as.character(df[[p]])
          v[is.na(v) | v == "" | v == "NA"] <- NA
          v
        })
        mat <- do.call(cbind, parts)
        split_group <- apply(mat, 1, function(r) {
          r <- r[!is.na(r)]
          if (length(r) == 0) "n/a" else paste(r, collapse = " | ")
        })
        df$color_group <- factor(split_group, levels = sort(unique(split_group)))
        color_by <- paste(present, collapse = " | ")
      } else {
        df$color_group <- factor(df$dataset, levels = sort(unique(df$dataset)))
        color_by <- "dataset"
      }
      attr(df, "color_by") <- color_by
      attr(df, "has_split") <- has_split
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

