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
                                          reload_pgxdir,
                                          current_page) {
  shiny::moduleServer(
    id, function(input, output, session) {
      ns <- session$ns ## NAMESPACE

      nr_ds_received <- reactiveVal(0)

      # callbackR for the "New dataset received" modal: jump to Shared datasets
      show_shared_tab <- function(value) {
        if (isTRUE(value)) {
          bigdash.selectTab(session, "sharing-tab")
        }
      }

      # keep navbar badges in sync with the count of pending received datasets
      update_shared_badges <- function(n) {
        shinyjs::runjs(sprintf(
          "(function(n){
            var sub = document.querySelector('a[data-target=\"sharing-tab\"]');
            if (sub) {
              var oldSub = sub.querySelector('.shared-pending-badge');
              if (oldSub) oldSub.remove();
              if (n > 0) {
                var b = document.createElement('span');
                b.className = 'shared-pending-badge';
                b.style.cssText = 'display:inline-block;background:#dc3545;color:#fff;border-radius:10px;padding:0 6px;margin-left:6px;font-size:11px;font-weight:600;line-height:16px;vertical-align:middle;';
                b.textContent = n;
                sub.appendChild(b);
              }
            }
            document.querySelectorAll('.nav-link.dropdown-toggle').forEach(function(el){
              if (el.textContent.trim().indexOf('Datasets') !== 0) return;
              var oldDot = el.querySelector('.shared-pending-dot');
              if (oldDot) oldDot.remove();
              if (n > 0) {
                var d = document.createElement('span');
                d.className = 'shared-pending-dot';
                d.style.cssText = 'display:inline-block;background:#dc3545;border-radius:50%%;width:8px;height:8px;margin-left:6px;vertical-align:middle;';
                el.appendChild(d);
              }
            });
          })(%d);", n
        ))
      }

      ## ------------ get received files
      getReceivedFiles <- shiny::reactivePoll(
        intervalMillis = 10000,
        session = NULL,
        checkFunc = function() {
          if (!auth$logged || auth$email == "") {
            return(nr_ds_received())
          }
          current_user <- auth$email
          pgxfiles <- dir(
            path = pgx_shared_dir,
            pattern = paste0("__to__", current_user, "__from__.*__$"),
            ignore.case = TRUE
          )

          current_ds_received <- length(pgxfiles)
          if (length(pgxfiles) > nr_ds_received()) {
            # modal that tells that user received a new dataset
            shinyalert::shinyalert(
              "New dataset received!",
              paste(
                "You have received a dataset from another user.",
                "Click below to view and accept it."
              ),
              showConfirmButton = TRUE,
              confirmButtonText = "Go to shared datasets",
              confirmButtonCol = "#337ab7",
              showCancelButton = TRUE,
              cancelButtonText = "Later",
              callbackR = show_shared_tab
            )
          }
          nr_ds_received(current_ds_received)
          return(nr_ds_received())
        },
        valueFunc = function() {
          req(auth$logged)
          if (!auth$logged || auth$email == "") {
            return(NULL)
          }

          # get received pgx files
          current_user <- auth$email
          pgxfiles <- dir(
            path = pgx_shared_dir,
            pattern = paste0("__to__", current_user, "__from__.*__$"),
            ignore.case = TRUE
          )
          return(pgxfiles)
        }
      )

      # refresh navbar badges whenever the pending-shares count changes
      shiny::observe({
        files <- getReceivedFiles()
        update_shared_badges(if (is.null(files)) 0L else length(files))
      })

      receivedPGXtable <- shiny::eventReactive(
        c(getReceivedFiles()),
        {
          received_files <- getReceivedFiles()
          if (is.null(received_files) || length(received_files) == 0) {
            df <- data.frame(
              Dataset = "-",
              To = "-",
              Actions = "-"
            )
            ## return(NULL)
          } else {
            # split the file name into user who shared and file name
            received_pgx <- sub("__to__.*", "", received_files)
            received_from <- gsub(".*__from__|__$", "", received_files)

            accept_btns <- makebuttonInputs2(
              FUN = actionButton,
              len = received_files,
              id = ns("accept_pgx__"),
              label = "",
              width = "50px",
              inline = TRUE,
              icon = shiny::icon("check"),
              class = "btn-inline btn-success",
              style = "padding:0px; margin:0px; font-size:85%;",
              tooltip = "Accept dataset",
              ## onclick = paste0('Shiny.onInputChange(\"',ns("accept_pgx"),'\", this.id, {priority: "event"})')
              onclick = paste0('Shiny.onInputChange("', ns("accept_pgx"), '", this.id, {priority: "event"})')
            )

            decline_btns <- makebuttonInputs2(
              FUN = actionButton,
              len = received_files,
              id = "decline_pgx__",
              label = "",
              width = "50px",
              inline = TRUE,
              icon = shiny::icon("x"),
              class = "btn-inline btn-danger",
              style = "padding:0px; margin:0px; font-size:85%;",
              tooltip = "Declinie dataset",
              onclick = paste0('Shiny.onInputChange(\"', ns("decline_pgx"), '\", this.id, {priority: "event"})')
            )

            df <- data.frame(
              Dataset = received_pgx,
              From = received_from,
              Actions = paste(accept_btns, decline_btns)
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
            DT::formatStyle(0, target = "row", fontSize = "14px", lineHeight = "90%")
        }
      )

      # ----------------- event when a shared pgx is accepted by a user
      observeEvent(input$accept_pgx,
        {
          # get pgx name and remove the __from__* tag
          pgx_name <- stringr::str_split(input$accept_pgx, "accept_pgx__")[[1]][2]
          new_pgx_name <- stringr::str_split(pgx_name, "__from__")[[1]][1]
          new_pgx_name <- sub("__to__.*", "", pgx_name)

          ## rename the file to be a valid pgx file
          pgxdir <- auth$user_dir
          file_from <- file.path(pgx_shared_dir, pgx_name)
          file_to <- file.path(pgxdir, new_pgx_name)

          ## check number of datasets
          numpgx <- length(dir(pgxdir, pattern = "*.pgx$"))
          if (!auth$options$ENABLE_DELETE) numpgx <- length(dir(pgxdir, pattern = "*.pgx$|*.pgx_$"))
          maxpgx <- as.integer(auth$options$MAX_DATASETS)
          if (numpgx >= maxpgx) {
            shinyalert_storage_full(numpgx, maxpgx, auth$level) ## from ui-alerts.R
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
          if (!file.rename(file_from, file_to)) {
            info("[loading_server.R] accept_pgx : rename failed. trying file.copy ")
            ## file.rename does not allow "cross-device link"
            file.copy(file_from, file_to)
            file.remove(file_from)
          }

          ## reload pgx dir so the newly accepted pgx files are registered in user table
          reload_pgxdir(reload_pgxdir() + 1)
        },
        ignoreInit = TRUE
      )

      #--------------  event when a shared pgx is declined by a user
      observeEvent(input$decline_pgx,
        {
          pgx_name <- stringr::str_split(input$decline_pgx, "decline_pgx__")[[1]][2]
          received_file <- file.path(pgx_shared_dir, pgx_name)
          file.remove(received_file)
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
