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

    ns <- session$ns

    pgxTable_DT <- reactive({
      df <- rl$pgxTable_data

      # need this, otherwise there is an error on user logout
      if (length(df$dataset) == 0) df <- NULL

      req(df)

      target1 <- grep("date", colnames(df))
      target2 <- grep("description", colnames(df))
      target3 <- grep("conditions", colnames(df))
      target4 <- grep("dataset", colnames(df))

      # create action menu for each row
      menus <- c()
      for (i in 1:nrow(df)) {
        new_menu <- actionMenu(
          div(
            style = "width: 175px;",
            div(
              shiny::actionButton(
                ns("deletebutton"),
                label = "Delete Dataset", icon = icon("trash"),
                class = "btn btn-outline-danger",
                width = '100%',
                style = 'border: none;'
              ),
              shiny::downloadButton(
                ns("downloadpgx"),
                label = "Download PGX",
                class = "btn btn-outline-dark",
                style = "width: 100%; border: none;"
              ),
              downloadButton2(
                ns("downloadzip"),
                label = "Download ZIP", icon = icon("file-archive"),
                class = "btn btn-outline-dark",
                style = "width: 100%; border: none;"
              ),
              shiny::actionButton(
                ns("sharebutton"),
                label = "Share Dataset", icon = icon("share-nodes"),
                class = "btn btn-outline-info",
                width = '100%',
                style = 'border: none;'
              ),
            )
          ),
          size = "sm",
          icon = shiny::icon("ellipsis-vertical"),
          status = "dark"
        )
        menus <- c(menus, as.character(new_menu))
      }

      df$actions <- menus
      colnames(df)[ncol(df)] <- ' '

      DT::datatable(
        df,
        class = "compact hover",
        rownames = TRUE,
        escape = FALSE,
        editable = list(
          target = 'cell',
          disable = list(columns = c(1,3:ncol(df)))
        ),
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
            list(width = "30vw", targets = target2),
            list(sortable = FALSE, targets = ncol(df))
          )
        ) ## end of options.list
      )
    })

    # make changes to pgxtable
    observeEvent(
      input[['datasets-datatable_cell_edit']], {
        row <- input[['datasets-datatable_cell_edit']]$row
        col <- input[['datasets-datatable_cell_edit']]$col
        val <- input[['datasets-datatable_cell_edit']]$value
        rl$pgxTable_data[row, col] <- val
        rl$pgxTable_edited <- rl$pgxTable_edited + 1
        rl$pgxTable_edited_row <- row
        rl$pgxTable_edited_col <- col
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
