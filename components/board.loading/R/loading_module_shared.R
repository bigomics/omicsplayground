##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

upload_module_shared_ui <- function(id, height = 720) {
  ns <- shiny::NS(id)
}

upload_module_shared_server <- function(id,
                                        auth,
                                        pgx_shared_dir,
                                        sendShareMessage,
                                        current_page,
                                        refresh) {
  shiny::moduleServer(
    id, function(input, output, session) {
      ns <- session$ns ## NAMESPACE

      refresh_table <- reactiveVal(0)

      ## ------------ get received files
      getSharedFiles <- shiny::reactive({
        req(auth$logged)
        if (!auth$logged) {
          return(c())
        }
        if (auth$email == "") {
          return(c())
        }
        ## allow trigger for when a shared pgx is accepted / decline
        refresh_table()
        refresh()
        pgxfiles <- dir(
          path = pgx_shared_dir,
          pattern = paste0("__to__.*__from__", auth$email, "__$"),
          ignore.case = TRUE
        )
        return(pgxfiles)
      })

      sharedPGXtable <- shiny::eventReactive(
        c(current_page(), getSharedFiles()),
        {
          ## req(current_page() == "sharing-tab")
          shared_files <- getSharedFiles()
          if (length(shared_files) == 0) {
            df <- data.frame(
              Dataset = "-",
              To = "-",
              Actions = "-"
            )
            ## return(NULL)
          } else {
            # split the file name into user who shared and file name
            shared_pgx <- sub("__to__.*", "", shared_files)
            shared_to <- gsub(".*__to__|__from__.*", "", shared_files)
            shared_from <- gsub(".*__from__|__$", "", shared_files)

            resend_btns <- makebuttonInputs2(
              FUN = actionButton,
              len = shared_files,
              id = ns("resend_pgx__"),
              label = "",
              width = "50px",
              inline = TRUE,
              icon = shiny::icon("repeat"),
              class = "btn-inline btn-success",
              style = "padding:0px; margin:0px; font-size:85%;",
              tooltip = "Resend this dataset",
              onclick = paste0('Shiny.onInputChange("', ns("resend_pgx"), '", this.id, {priority: "event"})')
            )

            cancel_btns <- makebuttonInputs2(
              FUN = actionButton,
              len = shared_files,
              id = "cancel_pgx__",
              label = "",
              width = "50px",
              inline = TRUE,
              icon = shiny::icon("x"),
              class = "btn-inline btn-danger",
              style = "padding:0px; margin:0px; font-size:85%;",
              tooltip = "Cancel this dataset sharing",
              onclick = paste0('Shiny.onInputChange(\"', ns("cancel_pgx"), '\", this.id, {priority: "event"})')
            )

            df <- data.frame(
              Dataset = shared_pgx,
              To = shared_to,
              Actions = paste(resend_btns, cancel_btns)
            )
          }

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
            DT::formatStyle(0, target = "row", fontSize = "14px", lineHeight = "95%")
        }
      )

      # ----------------- event when a shared pgx is accepted by a user
      observeEvent(input$resend_pgx,
        {
          pgx_name <- stringr::str_split(input$resend_pgx, "resend_pgx__")[[1]][2]
          sender <- gsub(".*__from__|__$", "", pgx_name)
          share_user <- gsub(".*__to__|__from__.*", "", pgx_name)
          gmail_creds <- file.path(ETC, "gmail_creds")
          sendShareMessage(pgx_name, sender, share_user, path_to_creds = gmail_creds)
        },
        ignoreInit = TRUE
      )

      #--------------  event when a shared pgx is declined by a user
      observeEvent(input$cancel_pgx,
        {
          pgx_name <- stringr::str_split(input$cancel_pgx, "cancel_pgx__")[[1]][2]
          shared_file <- file.path(pgx_shared_dir, pgx_name)
          file.remove(shared_file)
          refresh_table(refresh_table() + 1)
        },
        ignoreInit = TRUE
      )


      ## list of reactive objects acts like API or public function interface
      rlist <- list(
        getSharedFiles = getSharedFiles,
        sharedPGXtable = sharedPGXtable
      )

      return(rlist)
    }
  ) ## end of moduleServer
}
