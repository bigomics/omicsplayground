##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

upload_module_received_ui <- function(id, height = 720) {
  ns <- shiny::NS(id)
}


upload_module_received_server <- function(id,
                                          auth,
                                          pgx_shared_dir,
                                          getPGXDIR,
                                          max_datasets,
                                          enable_delete = TRUE,
                                          r_global) {
  shiny::moduleServer(
    id, function(input, output, session) {
      ns <- session$ns ## NAMESPACE

      refresh_table <- reactiveVal(0)

      ## ------------ get received files
      getReceivedFiles <- shiny::reactive({
        req(auth)
        if (!auth$logged()) {
          return(c())
        }
        if (auth$email() == "") {
          return(c())
        }
        ## allow trigger for when a shared pgx is accepted / decline
        refresh_table()
        pgxfiles <- dir(
          path = pgx_shared_dir,
          pattern = paste0("__to__", auth$email(), "__from__.*__$")
        )
        return(pgxfiles)
      })

      makebuttonInputs2 <- function(FUN, len, id, ...) {
        inputs <- character(length(len))
        for (i in seq_along(len)) {
          inputs[i] <- as.character(FUN(paste0(id, len[i]), ...))
        }
        inputs
      }

      receivedPGXtable <- shiny::eventReactive(
        c(r_global$nav, getReceivedFiles()),
        {
          req(r_global$nav == "load-tab")
          shared_files <- getReceivedFiles()
          if (length(shared_files) == 0) {
            return(NULL)
          }

          # split the file name into user who shared and file name
          shared_pgx <- sub("__to__.*", "", shared_files)
          shared_from <- gsub(".*__from__|__$", "", shared_files)

          accept_btns <- makebuttonInputs2(
            FUN = actionButton,
            len = shared_files,
            id = ns("accept_pgx__"),
            label = "",
            width = "50px",
            inline = TRUE,
            icon = shiny::icon("check"),
            class = "btn-inline btn-success",
            style = "padding:0px; margin:0px; font-size:85%;",

            ## onclick = paste0('Shiny.onInputChange(\"',ns("accept_pgx"),'\", this.id, {priority: "event"})')
            onclick = paste0('Shiny.onInputChange("', ns("accept_pgx"), '", this.id, {priority: "event"})')
          )

          decline_btns <- makebuttonInputs2(
            FUN = actionButton,
            len = shared_files,
            id = "decline_pgx__",
            label = "",
            width = "50px",
            inline = TRUE,
            icon = shiny::icon("x"),
            class = "btn-inline btn-danger",
            style = "padding:0px; margin:0px; font-size:85%;",
            onclick = paste0('Shiny.onInputChange(\"', ns("decline_pgx"), '\", this.id, {priority: "event"})')
          )

          df <- data.frame(
            Dataset = shared_pgx,
            From = shared_from,
            Actions = paste(accept_btns, decline_btns)
          )

          dt_table <- DT::datatable(
            df,
            rownames = FALSE,
            escape = FALSE,
            selection = "none",
            class = "compact row-border",
            options = list(
              dom = "t",
              pageLength = 999
            )
          ) %>%
            DT::formatStyle(0, target = "row", fontSize = "14px", lineHeight = "90%")
        }
      )

      # ----------------- event when a shared pgx is accepted by a user
      observeEvent(input$accept_pgx,
        {
          dbg("[loading_module_usershare:observeEvent(input$accept_pgx)] reacted!")

          # get pgx name and remove the __from__* tag
          pgx_name <- stringr::str_split(input$accept_pgx, "accept_pgx__")[[1]][2]
          new_pgx_name <- stringr::str_split(pgx_name, "__from__")[[1]][1]
          new_pgx_name <- sub("__to__.*", "", pgx_name)

          ## rename the file to be a valid pgx file
          pgxdir <- getPGXDIR()
          file_from <- file.path(pgx_shared_dir, pgx_name)
          file_to <- file.path(pgxdir, new_pgx_name)

          ## check number of datasets
          numpgx <- length(dir(pgxdir, pattern = "*.pgx$"))
          if (!enable_delete) numpgx <- length(dir(pgxdir, pattern = "*.pgx$|*.pgx_$"))
          maxpgx <- as.integer(max_datasets)
          if (numpgx >= maxpgx) {
            ## should use sprintf or glue here...

            msg <- "You have reached your datasets limit. Please delete some datasets, or <a href='https://events.bigomics.ch/upgrade' target='_blank'><b><u>UPGRADE</u></b></a> your account."

            shinyalert::shinyalert(
              title = "Your storage is full",
              text = HTML(msg),
              html = TRUE,
              type = "warning"
            )
            return(NULL)
          }

          ## warning on overwrite
          if (file.exists(file_to)) {
            shinyalert::shinyalert(
              "File already exists!",
              paste(
                "You have already a dataset called '", new_pgx_name,
                "'. Please delete it before accepting the new file."
              ),
              confirmButtonText = "Cancel",
              showCancelButton = FALSE
            )
            return(NULL)
          }

          ## Rename file to user folder. Some servers do not allow
          ## "cross-device link" and then we resort to slower copy
          dbg("[loading_server.R] accept_pgx : renaming file from = ", file_from, "to = ", file_to)
          if (!file.rename(file_from, file_to)) {
            info("[loading_server.R] accept_pgx : rename failed. trying file.copy ")
            ## file.rename does not allow "cross-device link"
            file.copy(file_from, file_to)
            file.remove(file_from)
          }

          ## reload pgx dir so the newly accepted pgx files are registered in user table
          r_global$reload_pgxdir <- r_global$reload_pgxdir + 1

          ## remove the accepted pgx from the table
          refresh_table(refresh_table() + 1)
        },
        ignoreInit = TRUE
      )

      #--------------  event when a shared pgx is declined by a user
      observeEvent(input$decline_pgx,
        {
          pgx_name <- stringr::str_split(input$decline_pgx, "decline_pgx__")[[1]][2]
          shared_file <- file.path(pgx_shared_dir, pgx_name)
          file.remove(shared_file)

          # remove the declined pgx from the table
          refresh_table(refresh_table() + 1)
        },
        ignoreInit = TRUE
      )


      ## list of reactive objects acts like API or public function interface
      rlist <- list(
        getReceivedFiles = getReceivedFiles,
        receivedPGXtable = receivedPGXtable
      )

      return(rlist)
    }
  ) ## end of moduleServer
}
