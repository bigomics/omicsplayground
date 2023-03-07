loading_table_datasets_ui <- function(id, height, width) {
  ns <- shiny::NS(id)

  info_text <- "This table contains a general information about all available datasets within the platform. For each dataset, it reports a brief description as well as the total number of samples, genes, gene sets (or pathways), corresponding phenotypes and the creation date."

  TableModuleUI(
    ns("datasets"),
    info.text = info_text,
    width = width,
    height = height,
    title = "Data files"
  )
}

loading_table_datasets_server <- function(id,
                                          rl) {
  moduleServer(id, function(input, output, session) {

    pgxTable_DT <- function() {
      df <- rl$pgxTable_data

      # need this, otherwise there is an error on user logout
      if (length(df$dataset) == 0) df <- NULL

      req(df)

      target1 <- grep("date", colnames(df))
      target2 <- grep("description", colnames(df))
      target3 <- grep("conditions", colnames(df))
      target4 <- grep("dataset", colnames(df))

      DT::datatable(
        df,
        class = "compact hover",
        rownames = TRUE,
        editable = TRUE,
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        fillContainer = TRUE,
        options = list(
          dom = "ft",
          pageLength = 9999,
          scrollX = FALSE,
          scrollY = "55vh",
          deferRender = TRUE,
          autoWidth = TRUE,
          columnDefs = list(
            list(width = "60px", targets = target1),
            list(width = "30vw", targets = target2)
          )
        ) ## end of options.list
      )
    }

    # make changes
    observeEvent(
      input[['datasets-datatable_cell_edit']], {
        row <- input[['datasets-datatable_cell_edit']]$row
        col <- input[['datasets-datatable_cell_edit']]$col
        val <- input[['datasets-datatable_cell_edit']]$value
        rl$pgxTable_data[row, col] <- val
        rl$pgxTable_edited <- rl$pgxTable_edited + 1
      }
    )

    pgxTable.RENDER <- function() {
      pgxTable_DT() %>%
        DT::formatStyle(0, target = "row", fontSize = "12px", lineHeight = "95%")
    }

    pgxTable_modal.RENDER <- function() {
      pgxTable_DT() %>%
        DT::formatStyle(0, target = "row", fontSize = "16px", lineHeight = "95%")
    }

    TableModuleServer(
      "datasets",
      func = pgxTable.RENDER,
      func2 = pgxTable_modal.RENDER,
      selector = "single"
    )
  })
}
