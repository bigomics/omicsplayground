##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

loading_table_datasets_ui <- function(
  id,
  title,
  info.text,
  caption,
  height,
  width) {
  ns <- shiny::NS(id)


  options <- tagList(
      shiny::checkboxGroupInput(ns("flt_datatype"), "Datatype", choices = ""),
      shiny::checkboxGroupInput(ns("flt_organism"), "Organism", choices = "")
  )

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    caption = caption,
    width = width,
    height = height,
    title = title,
    options = options
  )
}

loading_table_datasets_server <- function(id,
                                          getPGXINFO,
                                          getPGXDIR,
                                          auth,
                                          rl,
                                          r_global,
                                          enable_pgxdownload = FALSE,
                                          enable_delete = FALSE,
                                          enable_public_share = TRUE,
                                          enable_user_share = TRUE
                                          ) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    getFilteredPGXINFO <- shiny::reactive({
        ## get the filtered table of pgx datasets
        req(auth)
        if (!auth$logged()) {
            warning("[LoadingBoard:getFilteredPGXINFO] user not logged in!
                    not showing table!")
            return(NULL)
        }
        df <- getPGXINFO()
        if (is.null(df)) {
            return(NULL)
        }

        pgxdir <- getPGXDIR()
        pgxfiles <- dir(pgxdir, pattern = ".pgx$")
        sel <- sub("[.]pgx$", "", df$dataset) %in% sub("[.]pgx$", "", pgxfiles)
        df <- df[sel, , drop = FALSE]

        ## Apply filters
        if (nrow(df) > 0) {
            f1 <- f2 <- f3 <- rep(TRUE, nrow(df))
            notnull <- function(x) !is.null(x) && length(x) > 0 && x[1] != "" && !is.na(x[1])
            if (notnull(input$flt_datatype)) f2 <- (df$datatype %in% input$flt_datatype)
            if (notnull(input$flt_organism)) f3 <- (df$organism %in% input$flt_organism)
            df <- df[which(f1 & f2 & f3), , drop = FALSE]
            df$date <- as.Date(df$date, format = "%Y-%m-%d")
            df <- df[order(df$date, decreasing = TRUE), ]
            if (nrow(df) > 0) rownames(df) <- nrow(df):1
        }

        kk <- unique(c(
            "dataset", "description", "datatype", "nsamples",
            "ngenes", "nsets", "conditions", "date", "organism",
            "creator"
        ))
        kk <- intersect(kk, colnames(df))
        df <- df[, kk, drop = FALSE]
        df
    })

    table_data <- shiny::reactive({
        #r_global$reload_pgxdir

        df <- getFilteredPGXINFO()

        df$dataset <- sub("[.]pgx$", "", df$dataset)
        df$conditions <- gsub("[,]", " ", df$conditions)
        df$conditions <- sapply(as.character(df$conditions), andothers, split = " ", n = 5)
        df$description <- playbase::shortstring(as.character(df$description), 200)
        df$nsets <- NULL
        df$organism <- NULL

        return(df)
    })

    andothers <- function(s, split = " ", n = 8) {
        if (is.na(s)) {
            return("")
        }
        s <- sub("^[ ]*", "", s)
        s <- sub("[ ]+", " ", s)
        s1 <- strsplit(s, split = split)[[1]]
        if (length(s1) <= n) {
            return(s)
        }
        n2 <- setdiff(length(s1), n)
        paste(paste(head(s1, n), collapse = " "), "(+", n2, "others)")
    }

    observeEvent(getPGXINFO(), {
        df <- getPGXINFO()
        if(is.null(df)) return()
        datatypes <- sort(setdiff(df$datatype, c(NA, "")))
        organisms <- sort(setdiff(df$organism, c(NA, "")))
        shiny::updateCheckboxGroupInput(session, "flt_datatype", choices = datatypes)
        shiny::updateCheckboxGroupInput(session, "flt_organism", choices = organisms)
    })

    pgxTable_DT <- reactive({

      df <- table_data()
      shiny::req(df)

      is.dt <- is.data.frame(df)
      if (!is.dt || nrow(df) == 0) {
        shinyalert::shinyalert(
          title = "Empty?",
          text = paste("Your dataset library seems empty. Please upload new data or import",
            "a dataset from the public datasets folder."
            )
        )
      }
      validate(need(nrow(df)>0, 'Need at least one dataset!'))

      ## need this, otherwise there is an error on user logout
      if (length(df$dataset) == 0) df <- NULL

      df$creator <- NULL
      target1 <- grep("date", colnames(df))
      target2 <- grep("description", colnames(df))
      target3 <- grep("conditions", colnames(df))
      target4 <- grep("dataset", colnames(df))

      # create action menu for each row
      menus <- c()
      for (i in 1:nrow(df)) {

        download_pgx_menuitem <- NULL
        share_public_menuitem <- NULL
        share_dataset_menuitem <- NULL
        delete_pgx_menuitem <- NULL

        if(enable_pgxdownload) {
          download_pgx_menuitem <- shiny::actionButton(
            ns(paste0("download_pgx_row_",i)),
            label = "Download PGX",
            icon = shiny::icon('download'),
            class = "btn btn-outline-dark",
            style = "border: none;",
            width = '100%',
            onclick=paste0('Shiny.onInputChange(\"',ns("download_pgx"),'\",this.id,{priority: "event"})')
          )
        }
        if(enable_public_share) {
          share_public_menuitem <- shiny::actionButton(
              ns(paste0("share_public_row_", i)),
              label = "Share Public",
              icon = shiny::icon('share-nodes'),
              class = "btn btn-outline-info",
              style = 'border: none;',
              width = '100%',
              onclick=paste0('Shiny.onInputChange(\"',ns("share_public_pgx"),'\",this.id,{priority: "event"})')
          )
        }
        if(enable_user_share) {
          share_dataset_menuitem <- shiny::actionButton(
            ns(paste0("share_dataset_row_", i)),
            label = "Share with User",
            icon = shiny::icon('share-nodes'),
            class = "btn btn-outline-info",
            style = 'border: none;',
            width = '100%',
            onclick=paste0('Shiny.onInputChange(\"',ns("share_pgx"),'\",this.id,{priority: "event"})')
          )
        }

        if(enable_delete) {
          delete_pgx_menuitem <- shiny::actionButton(
            ns(paste0("delete_dataset_row_",i)),
            label = "Delete Dataset",
            icon = shiny::icon("trash"),
            class = "btn btn-outline-danger",
            style = 'border: none;',
            width = '100%',
            onclick=paste0('Shiny.onInputChange(\"',ns("delete_pgx"),'\",this.id,{priority: "event"});')
          )
        }

        new_menu <- actionMenu(  ## ui-DrowDownMenu.R
          div(
            style = "width: 160px;",
            div(
              download_pgx_menuitem,
              shiny::actionButton(
                ns(paste0("download_zip_row_", i)),
                label = "Download ZIP",
                icon = shiny::icon("file-archive"),
                class = "btn btn-outline-dark",
                style = "border: none;",
                width = '100%',
                onclick=paste0('Shiny.onInputChange(\"',ns("download_zip"),'\",this.id,{priority: "event"})')
                ),
              share_public_menuitem,
              share_dataset_menuitem,
              delete_pgx_menuitem
            )
          ),
          size = "sm",
          icon = shiny::icon("ellipsis-vertical"),
          status = "dark"
        )
        menus <- c(menus, as.character(new_menu))
      }

      observeEvent(input$download_pgx, { rl$download_pgx <- input$download_pgx })
      observeEvent(input$download_zip, { rl$download_zip <- input$download_zip })
      observeEvent(input$share_pgx, { rl$share_pgx <- input$share_pgx }, ignoreInit = TRUE)
      observeEvent(input$share_public_pgx, { rl$share_public_pgx <- input$share_public_pgx }, ignoreInit = TRUE)
      observeEvent(input$delete_pgx, { rl$delete_pgx <- input$delete_pgx }, ignoreInit = TRUE)

      DT::datatable(
        df,
        class = "compact hover",
        rownames = menus,
        escape = FALSE,
        editable = list(
          target = 'cell',
          disable = list(columns = c(1,3:ncol(df)))
        ),
        extensions = c("Scroller"),
        plugins = 'scrollResize',
        selection = list(mode = "single", target = "row", selected = 1),
        fillContainer = TRUE,
        options = list(
          dom = "ft",
          pageLength = 9999,
          scrollX = FALSE,
          scrollY = "55vh",
          scrollResize = TRUE,
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
        DT::formatStyle(0, target = "row", fontSize = "20px", lineHeight = "95%")
    }

    TableModuleServer(
      "datasets",
      func = pgxTable.RENDER,
      func2 = pgxTable_modal.RENDER,
      selector = "single"
    )

    return(
        pgxtable_data = table_data
    )
  })
}
