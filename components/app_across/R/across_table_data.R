##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' @export
across_table_data_ui <- function(id,
                                 title,
                                 height,
                                 width) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("table"),
    title = title,
    label = "d",
    info.text = "Data table showing values for selected features across all samples.",
    height = height,
    width = width
  )
}

#' @export
across_table_data_server <- function(id,
                                     getPlotData,
                                     dataset_info = shiny::reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    table_data <- shiny::reactive({
      df <- getPlotData()
      shiny::validate(shiny::need(
        !is.null(df) && nrow(df) > 0,
        "No data available. Please select features and click 'Query'."
      ))

      ## The value column header reflects the chosen value type / scale (Count,
      ## Z-score or log2); computed upstream in getPlotData where those inputs live.
      value_label <- attr(df, "value_label")
      if (is.null(value_label)) value_label <- "Count"

      result <- data.frame(
        Feature = if ("gene" %in% colnames(df)) df$gene else NA,
        Dataset = if ("dataset" %in% colnames(df)) df$dataset else NA,
        stringsAsFactors = FALSE
      )

      ## Annotate each row with its dataset's datatype / organism when available
      ## (from datasets-info.csv via getDatasetInfoTileDB). Skipped if absent.
      info <- dataset_info()
      if (!is.null(info) && "dataset" %in% colnames(info)) {
        if ("datatype" %in% colnames(info)) {
          result$Datatype <- info$datatype[match(result$Dataset, info$dataset)]
        }
        if ("organism" %in% colnames(info)) {
          result$Organism <- info$organism[match(result$Dataset, info$dataset)]
        }
      }

      result$Sample <- if ("sample_short" %in% colnames(df)) df$sample_short else NA
      result[[value_label]] <- if ("count" %in% colnames(df)) round(df$count, 2) else NA

      ## "color_group" (and "split_group", if present) are derived grouping
      ## helpers; hide them. The selected phenotype columns still show.
      standard_cols <- c("gene", "dataset", "sample", "sample_short", "count", "color_group", "split_group")
      pheno_cols <- setdiff(colnames(df), standard_cols)
      for (col in pheno_cols) result[[col]] <- df[[col]]
      result$FullSampleID <- if ("sample" %in% colnames(df)) df$sample else NA

      attr(result, "value_col") <- value_label
      result
    })

    table.RENDER <- shiny::reactive({
      df <- table_data()
      value_col <- attr(df, "value_col")
      if (is.null(value_col)) value_col <- "Count"
      numeric_cols <- which(colnames(df) == value_col)

      DT::datatable(
        df,
        rownames = FALSE,
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row"),
        class = "compact cell-border stripe hover",
        plugins = "scrollResize",
        fillContainer = TRUE,
        options = list(
          dom = "frtip",
          pageLength = 100,
          lengthMenu = c(25, 40, 100, 250, 1000),
          scroller = TRUE,
          scrollY = SCROLLY_MODAL,
          scrollResize = TRUE,
          deferRender = TRUE
        )
      ) %>%
        DT::formatSignif(numeric_cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    })

    table.RENDER_modal <- shiny::reactive({
      dt <- table.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    TableModuleServer(
      "table",
      func = table.RENDER,
      func2 = table.RENDER_modal,
      csvFunc = table_data,
      selector = "none"
    )
  })
}

