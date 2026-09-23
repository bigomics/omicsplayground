##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

LoadingBoard <- function(id,
                         pgx,
                         auth,
                         limits = c(
                           "samples" = 1000, "comparisons" = 20,
                           "genes" = 20000, "genesets" = 10000,
                           "datasets" = 10
                         ),
                         pgx_topdir,
                         load_example,
                         load_example_dataset = NULL,
                         reload_pgxdir,
                         current_page,
                         load_uploaded_data,
                         recompute_pgx,
                         new_upload,
                         save_pgx = NULL,
                         pgx_save_target = NULL,
                         pgx_source_dir = NULL,
                         parent) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE


    reload_pgxdir_public <- reactiveVal(0)
    refresh_shared <- reactiveVal(0)
    is_data_loaded <- reactiveVal(NULL)

    ## static, not changing
    pgx_shared_dir <- stringr::str_replace_all(pgx_topdir,
                                               c("data" = "data_shared"))
    pgx_public_dir <- stringr::str_replace_all(pgx_topdir,
                                               c("data" = "data_public"))
    enable_public_tabpanel <- dir.exists(pgx_public_dir)

    ## -------------------------------------------------------------------
    ## Received/shared UI
    ## -------------------------------------------------------------------

    ## Navigate to the "Shared datasets" inner tab of the Library board. The old
    ## bigdash "sharing-tab" board was folded into this tabset during the UI
    ## reshuffle, so bigdash.selectTab() no longer reaches it.
    goto_sharing_tab <- function() {
      bslib::nav_select("app-sidebar", "Library", session = parent)
      shiny::updateTabsetPanel(session, "tabs", selected = "sharing_tab")
    }

    pgxreceived <- upload_module_received_server(
      id = "received",
      auth = auth,
      pgx_shared_dir = pgx_shared_dir,
      ##      max_datasets = auth$options$MAX_DATASETS,  ## wrong and not needed...
      reload_pgxdir = reload_pgxdir,
      current_page = current_page,
      goto_sharing_tab = goto_sharing_tab
    )

    pgxshared <- upload_module_shared_server(
      id = "shared",
      auth = auth,
      pgx_shared_dir = pgx_shared_dir,
      sendShareMessage = sendShareMessage,
      current_page = current_page,
      refresh = refresh_shared
    )

    output$sharing_alert <- renderUI({
      received_files <- pgxreceived$getReceivedFiles()
      shared_files <- pgxshared$getSharedFiles()
      num_received <- length(received_files)
      num_shared <- length(shared_files)

      no_sharing1 <- !auth$options$ENABLE_USER_SHARE
      no_sharing2 <- (num_received == 0 && num_shared == 0)
      no_sharing <- no_sharing1 || no_sharing2

      if (no_sharing) {
        tag <- bs_alert(HTML(
          "This table shows the <b>available datasets</b> in your library. ",
          "The <b>Signature t-SNE</b> shows similarity clustering of ",
          "signatures using t-SNE. Select a dataset in the table and click ",
          "the <b>Load selected</b> button below."
        ))
        return(tag)
      }

      ## If not show alerts for sharing
      msg <- c()
      if (num_received > 0) {
        msg <- paste("You have received <strong>", num_received,
                     "datasets</strong> that you need to accept.")
      }
      if (num_shared > 0) {
        msg1 <- paste("You have still <strong>", num_shared,
                      "shared datasets</strong> waiting in the queue.")
        msg <- c(msg, msg1)
      }
      bs_alert(
        style = "warning",
        conditional = FALSE,
        shiny::HTML(paste(msg, "Please check the Sharing panel."))
      )
    })

    ## ======================================================================
    ## LOAD EXAMPLE TRIGGER
    ## ======================================================================
    observeEvent(load_example(), {

      ## Which dataset counts as "the example" is caller-chosen (e.g. the
      ## MultiOmics dashboard's popup wants "mox-brca" instead of the
      ## default "example-data") -- isolate() since this reactiveVal is only
      ## meant to be read at the moment load_example() itself fires, not to
      ## add its own dependency here.
      example_name <- if (!is.null(load_example_dataset)) {
        shiny::isolate(load_example_dataset())
      } else {
        "example-data"
      }

      # get the row which corresponds to the target example dataset
      data_names <- as.character(pgxtable$data()$dataset)
      example_row <- which(data_names == example_name)[1]
      has.exampledata <- (example_name %in% data_names)

      # if not found, throw error modal that the example dataset doesnt exist
      ## if (is.na(example_row)) {
      if (!has.exampledata) {
        shinyalert::shinyalert(
          title = "No example data",
          text = paste0("Sorry, the example dataset '", example_name, "' could not be found."),
          type = "warning",
          closeOnClickOutside = FALSE
        )
        return(NULL)
      } else {
        loadAndActivatePGX(example_name)
      }
    },
    ignoreInit = TRUE
    )

    ## ================================================================================
    ## Modules
    ## ================================================================================

    pgxtable <- loading_table_datasets_server(
      id = "pgxtable",
      pgx_shared_dir = pgx_shared_dir,
      pgx_public_dir = pgx_public_dir,
      pgx_topdir = pgx_topdir,
      auth = auth,
      loadAndActivatePGX = loadAndActivatePGX,
      loadPGX = loadPGX,
      refresh_shared = refresh_shared,
      reload_pgxdir_public = reload_pgxdir_public,
      reload_pgxdir = reload_pgxdir,
      recompute_pgx = recompute_pgx,
      loadbutton = reactive(input$loadbutton),
      new_upload = new_upload
    )

    loading_tsne_server(
      id = "tsne",
      pgx.dirRT = reactive(auth$user_dir),
      info.table = reactive(pgxtable$data()),
      r_selected = reactive(pgxtable$rows_all()),
      watermark = WATERMARK
    )


    if (enable_public_tabpanel) {
      pgxtable_public <- loading_table_datasets_public_server(
        id = "pgxtable_public",
        pgx_public_dir = pgx_public_dir,
        reload_pgxdir_public = reload_pgxdir_public,
        auth = auth,
        reload_pgxdir = reload_pgxdir,
        loadAndActivatePGX = loadAndActivatePGX
      )

      loading_tsne_server(
        id = "tsne_public",
        pgx.dir = reactive(pgx_public_dir),
        info.table = reactive(pgxtable_public$data()),
        r_selected = reactive(pgxtable_public$rows_all()),
        watermark = WATERMARK
      )
    }

    shiny::observeEvent(auth, {
      pgx_archive_dir <- file.path(auth$user_dir, "data_archive")
      enable_archive_tabpanel <- dir.exists(pgx_archive_dir)

      ## Show/hide archive tab based on whether directory exists
      if (enable_archive_tabpanel) {
        shiny::showTab(inputId = "tabs", target = "archive_tab", session = session)
        pgxtable_archive <- loading_table_datasets_public_server(
          id = "pgxtable_archive",
          pgx_public_dir = pgx_archive_dir,
          reload_pgxdir_public = reload_pgxdir_public,
          auth = auth,
          reload_pgxdir = reload_pgxdir,
          loadAndActivatePGX = loadAndActivatePGX
        )

        loading_tsne_server(
          id = "tsne_archive",
          pgx.dir = reactive(pgx_archive_dir),
          info.table = reactive(pgxtable_archive$data()),
          r_selected = reactive(pgxtable_archive$rows_all()),
          watermark = WATERMARK
        )
      } else {
        shiny::hideTab(inputId = "tabs", target = "archive_tab", session = session)
      }
    })

    ## -----------------------------------------------------------------------------
    ## Description
    ## -----------------------------------------------------------------------------

    shiny::observeEvent(input$module_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Loading a dataset from your library</strong>"),
        shiny::HTML(module_infotext),
        easyClose = TRUE, size = "xl"
      ))
    })

    module_infotext <- tspan(paste0(
      "This panel shows the available datasets within the platform. The table
        reports a brief description as well as the total number of samples,
        genes, gene sets (or pathways), corresponding phenotypes and the creation
        date.<br><br><b>Selecting the dataset:</b> Users can select a dataset in
        the table. The Dataset info shows the information of the dataset of
        interest and users can analyze the data by clicking the 'Analyze dataset'
        button.<br><br><br><center><iframe width='560' height='315'
        src='https://www.youtube.com/embed/elwT6ztt3Fo'
        title='YouTube video player' frameborder='0'
        allow='accelerometer; autoplay; clipboard-write; encrypted-media;
        gyroscope; picture-in-picture' allowfullscreen></iframe><center>"
    ), js = FALSE)
    module_infotext <- paste0(
      "<center><iframe width='560' height='315' src='https://www.youtube.com/embed/YTzLkio4M_4?si=LljECgKnb0TsgZDZ&amp;start=517' title='YouTube video player' frameborder='0' allow='accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share' referrerpolicy='strict-origin-when-cross-origin' allowfullscreen></iframe></center>"
    )

    sharing_infotext <- tspan(paste0(
      "<center><iframe width='560' height='315' src='https://www.youtube.com/embed/YTzLkio4M_4?si=oF52aFP_1knRilod&amp;start=551' title='YouTube video player' frameborder='0' allow='accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share' referrerpolicy='strict-origin-when-cross-origin' allowfullscreen></iframe></center>"
    ))

    shiny::observeEvent(input$loading_sharing, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Sharing a dataset from your library</strong>"),
        shiny::HTML(sharing_infotext),
        easyClose = TRUE, size = "xl"
      ))
    })

    ## =============================================================================
    ## ========================== OBSERVE/REACT ====================================
    ## =============================================================================

    ## =========================== BUTTON ACTIONS =============================
    ## disable button if no row is selected
    observeEvent(pgxtable$rows_selected(),
      {
        shiny::req(pgxtable)
        if (is.null(pgxtable$rows_selected())) {
          shinyjs::disable(id = "loadbutton")
        } else {
          shinyjs::enable(id = "loadbutton")
        }
      },
      ignoreNULL = FALSE
    )

    loadPGX <- function(pgxfile, pgxdir = NULL) {
      req(auth$logged)
      if (!auth$logged) {
        return(NULL)
      }

      ## Use provided directory or default to user directory
      if (is.null(pgxdir)) {
        pgxdir <- auth$user_dir
      }

      pgxfile <- paste0(sub("[.]pgx$", "", pgxfile), ".pgx") ## add/replace .pgx
      pgxfile1 <- file.path(pgxdir, pgxfile)

      pgx <- NULL
      if (file.exists(pgxfile1)) {
        pgx <- playbase::pgx.load(pgxfile1)
      } else {
        warning("[LoadingBoard::loadPGX] ***ERROR*** file not found : ", pgxfile1)
        return(NULL)
      }
      if (!is.null(pgx)) {
        pgx$name <- pgxfile
        return(pgx)
      } else {
        warning("[LoadingBoard::loadPGX] ***ERROR*** loading pgx object")
        return(NULL)
      }
    }

    ## Stable per-user half of the run key. Two tabs of the same user on the
    ## same dataset must join one run; two users must not. The save path
    ## already encodes whose copy is being written, so an anonymous deployment
    ## can safely collapse onto one literal - AI Studio falls back to the same
    ## one (app_studio/R/aireport_server.R), and both halves must agree or the
    ## two paths would key the same run differently.
    ai_report_user_key <- function() {
      email <- auth$email
      if (!is.null(email) && length(email) == 1L && !is.na(email) &&
            nzchar(email)) {
        return(as.character(email))
      }
      "anonymous"
    }

    ## Progress bar and completion handling for one subscription to a report
    ## run. Shared by the "start a run" and the "join the run already in
    ## flight" paths, so a session that joins gets the same feedback as the
    ## one that started it.
    ##
    ## Everything that talks back to the browser is bound to this session's
    ## reactive domain: promises restores the domain of the session that
    ## REGISTERED the run inside these callbacks, so a joiner's alert would
    ## otherwise be resolved through getDefaultReactiveDomain() and delivered
    ## to the starter - which also never receives shinyalert's JS dependencies
    ## for the joining session, since shinyalert inserts them into the default
    ## domain regardless of its own `session` argument.
    ai_report_run_handlers <- function(pgx_list, token, status = NULL) {
      progress <- shiny::Progress$new(session, min = 0, max = 1)
      joined <- !is.null(status) && !is.na(status$total) && status$total > 0L
      progress$set(
        message = "Generating AI reports",
        value = if (joined) status$done / status$total else 0,
        detail = if (joined) {
          paste0("Joined run in progress (", status$done, "/", status$total, ")")
        } else {
          "Preparing prompts..."
        }
      )
      closed <- FALSE
      close_progress <- function() {
        if (closed) return(invisible(NULL))
        closed <<- TRUE
        tryCatch(progress$close(), error = function(e) NULL)
      }

      pgx_now <- function() shiny::isolate(shiny::reactiveValuesToList(pgx))
      ## The results belong to the pgx the run was started from. If the user
      ## has loaded something else meanwhile they must not be written over the
      ## dataset now on screen.
      still_current <- function() {
        identical(token, ai_report_dataset_token(pgx_now()))
      }

      show_progress <- function(done, total, slot, ok) {
        label <- ai_report_slot_label(slot, pgx_list)
        verb <- if (is.na(slot)) "Starting" else
          if (isTRUE(ok)) "Finished" else "Failed"
        tryCatch(progress$set(
          value = if (total > 0) done / total else 0,
          detail = paste0(verb, " ", label, " (", done, "/", total, ")")
        ), error = function(e) NULL)
      }

      show_done <- function(result) {
        close_progress()
        if (is.null(result$ai)) {
          shiny::showNotification("AI report generation failed.",
            type = "error", session = session)
          return(invisible(NULL))
        }
        ## The manager has already written the pgx on disk. Mirror the result
        ## into this session's reactive only when the dataset is still on
        ## screen, so the boards pick it up without a reload.
        if (!still_current()) {
          shiny::showNotification(
            "AI reports finished and were saved to their dataset.",
            type = "message", session = session)
          return(invisible(NULL))
        }
        shiny::isolate(ai_report_merge_into_reactive(pgx, result$ai))
        tryCatch(
          ai_telemetry_record_reports(pgx_now(), user_email = auth$email),
          error = function(e) NULL
        )
        shinyalert::shinyalert(
          title = "AI reports ready",
          text = if (result$failed > 0L) {
            paste0("Your AI reports are ready. ", result$failed,
              " section(s) could not be generated.")
          } else {
            "Your AI reports are ready."
          },
          type = if (result$failed > 0L) "warning" else "success",
          confirmButtonText = "OK",
          session = session
        )
        invisible(NULL)
      }

      list(
        close = close_progress,
        on_progress = function(done, total, slot, ok) {
          shiny::withReactiveDomain(session, show_progress(done, total, slot, ok))
        },
        on_done = function(result) {
          shiny::withReactiveDomain(session, show_done(result))
        }
      )
    }

    maybe_offer_ai_reports <- function(pgxfile, is_user_dir) {
      if (!isTRUE(opt$ENABLE_AI)) return(invisible(NULL))

      ## Only offer generation to users who can persist the result: dataset
      ## owners (loaded from their own dir) or, when the admin feature is
      ## enabled, admins (curators, who write back to the source dir via
      ## save_current_pgx). Otherwise the paid ai_report_generate() below runs
      ## only to no-op on save. Same gate as AI Studio on-demand generation.
      admin_ok <- isTRUE(auth$ADMIN) && isTRUE(opt$ENABLE_ADMIN)
      if (!isTRUE(is_user_dir) && !admin_ok) return(invisible(NULL))

      llm_model <- getUserOption(session, "llm_model")
      if (is.null(llm_model) || llm_model == "") return(invisible(NULL))

      pgx_list <- shiny::reactiveValuesToList(pgx)
      report_modules <- ai_report_modules_for_pgx(pgx_list)

      ## Token of the dataset the prompt refers to. If the user loads something
      ## else while generation is in flight, the results belong to a pgx that is
      ## no longer on screen and must not be written over the new one.
      token <- ai_report_dataset_token(pgx_list)
      ## Resolved here, while there is still session context: the run itself
      ## never touches the session again.
      save_path <- if (is.function(pgx_save_target)) {
        pgx_save_target(pgx_list)
      } else {
        NULL
      }
      run_key <- if (is.null(save_path)) NULL else {
        ai_report_run_key(save_path, token, ai_report_user_key(), llm_model)
      }

      ## A run already in flight owns this dataset's reports, so attach to it
      ## for progress rather than asking again. ai_report_needs_generation()
      ## reads the pgx, which does not change until the run commits at the very
      ## end - without this, every reload during a multi-minute run re-opens
      ## the prompt and invites a duplicate answer.
      if (!is.null(run_key) && ai_report_run_active(run_key)) {
        handlers <- ai_report_run_handlers(pgx_list, token,
          ai_report_run_status(run_key))
        sub_id <- ai_report_run_subscribe(run_key, handlers$on_progress,
          handlers$on_done)
        if (is.null(sub_id)) {
          ## Finished between the check and the subscribe.
          handlers$close()
          return(invisible(NULL))
        }
        session$onSessionEnded(function() {
          ai_report_run_unsubscribe(run_key, sub_id)
        })
        return(invisible(NULL))
      }

      ## Do not require every possible module report. Some modules are optional
      ## or can fail independently; any valid pgx$ai report is enough to avoid
      ## prompting on every load.
      if (!ai_report_needs_generation(pgx_list)) {
        return(invisible(NULL))
      }

      cred_fn <- get_ai_credentials(session)
      ds_name <- if (!is.null(pgx_list$name)) pgx_list$name else pgxfile

      shinyalert::shinyalert(
        title = "Missing AI reports",
        text = paste0("Dataset '", ds_name,
          "' has missing AI reports. Would you like to compute them now?",
          " This runs in the background - you can keep using the app."),
        type = "info",
        showCancelButton = TRUE,
        confirmButtonText = "Yes",
        cancelButtonText = "No",
        callbackR = function(confirmed) {
          if (!isTRUE(confirmed)) return(NULL)

          if (is.null(save_path)) {
            shiny::showNotification(
              "Cannot determine where to save this dataset; reports not generated.",
              type = "error", session = session)
            return(NULL)
          }

          ## The run is owned by the app, not by this session: closing the
          ## tab or loading another dataset no longer abandons work that has
          ## already been paid for, and a second tab on the same dataset joins
          ## the run in progress instead of starting a duplicate one. This
          ## session only subscribes for progress.
          handlers <- ai_report_run_handlers(pgx_list, token)
          run <- ai_report_run_start(
            pgx_list,
            save_path   = save_path,
            llm_model   = llm_model,
            select      = report_modules,
            credentials = cred_fn,
            user_key    = ai_report_user_key(),
            user_email  = auth$email,
            on_progress = handlers$on_progress,
            on_done     = handlers$on_done
          )

          if (is.null(run)) {
            handlers$close()
            return(NULL)
          }
          ## Unsubscribe by the id the manager handed back, not by the session
          ## token: answering "Yes" twice in one tab subscribes twice, and each
          ## subscription owns the progress bar only its own on_done closes.
          session$onSessionEnded(function() {
            ai_report_run_unsubscribe(run$key, run$sub_id)
          })
        }
      )
    }

    loadAndActivatePGX <- function(pgxfile, pgxdir = NULL) {

      ## During loading show loading pop-up modal
      firstpgx <- (length(names(pgx))==0)
      if(firstpgx) {
        ui.showStartupModal()
      } else {
        ui.showCartoonModal()
      }

      loaded_pgx <- loadPGX(pgxfile, pgxdir = pgxdir)
      if (is.null(loaded_pgx)) {
        warning("[loadAndActivatePGX] ERROR loading PGX file ", pgxfile, "\n")
        beepr::beep(10)
        shiny::removeModal()
        return(NULL)
      }

      ## Record the source directory of the dataset just loaded so downstream
      ## save/authorization logic (e.g. admin write-back to public/shared
      ## datasets) knows where it actually lives.
      if (!is.null(pgx_source_dir)) {
        pgx_source_dir(if (is.null(pgxdir)) auth$user_dir else pgxdir)
      }

      ## ----------------- update PGX object ---------------------------------
      kk <- grep("name|date",names(loaded_pgx),invert=TRUE)
      size0 <- object.size(loaded_pgx[kk])
      shiny::withProgress(message = "Initializing. Please wait...", value = 0.33, {
        loaded_pgx <- playbase::pgx.initialize(loaded_pgx)

        if (is.null(loaded_pgx)) {
          warning("[loadAndActivatePGX] ERROR in object initialization\n")
          beepr::beep(10)
          shiny::showNotification("ERROR in object initialization!\n")
          shiny::removeModal()
          return(NULL)
        }
        loaded_pgx$name <- sub("[.]pgx$", "", pgxfile) ## always use filename

        ## if PGX object has been updated with pgx.initialize, we save
        ## the updated object (but only if loading from user directory)
        kk <- grep("name|date",names(loaded_pgx),invert=TRUE)
        size1 <- object.size(loaded_pgx[kk])
        is_user_dir <- is.null(pgxdir) || (pgxdir == auth$user_dir)
        if (size1 != size0 && is_user_dir && !is.null(save_pgx)) {
          info("[loadAndActivatePGX] WARNING: initialized PGX changed! saving updated PGX")
          save_pgx(loaded_pgx)
        }

        ## Copying to pgx list to reactiveValues in
        ## session environment.
        info("[loadAndActivatePGX] copying pgx object to global environment")
        empty.slots <- setdiff(names(pgx), names(loaded_pgx))
        isolate({
          for (e in empty.slots) {
            pgx[[e]] <- NULL
          }
          for (i in 1:length(loaded_pgx)) {
            pgx[[names(loaded_pgx)[i]]] <- loaded_pgx[[i]]
          }
        })
      }) ## end of withProgress

      ## clean up
      gc()
      remove(loaded_pgx)

      ## ----------------- AI reports: offer to compute if missing -----------
      maybe_offer_ai_reports(pgxfile, is_user_dir)

      ## notify new data uploaded
      if (is.null(is_data_loaded())) {
        is_data_loaded(1)
      } else {
        is_data_loaded(is_data_loaded() + 1)
      }

      info("[loadAndActivatePGX] done!")
    }

    observeEvent(input$newuploadbutton, {
      ##new_upload(new_upload() + 1)
      bslib::nav_select("app-sidebar", "Upload", session=parent)
    })

    observeEvent(load_uploaded_data(), {
      upload_pgx <- sub("[.]pgx$", "", load_uploaded_data())
      loadAndActivatePGX(upload_pgx)
      load_uploaded_data(NULL)
    })

    ## ================================================================================
    ## Header
    ## ================================================================================

    pgx_stats <- reactive({
      pgx_info <- pgxtable$data()
      shiny::req(pgx_info)
      ndatasets <- nrow(pgx_info)
      nsamples <- sum(as.integer(pgx_info$nsamples), na.rm = TRUE)
      FC.file <- file.path(auth$user_dir, "datasets-allFC.csv")
      if (file.exists(FC.file)) {
        contrasts <- get_contrasts_from_user(auth)
        ncontrasts <- sum(contrasts, na.rm = TRUE)
        return(
          paste(ndatasets, "Data sets &nbsp;&nbsp;&nbsp;", nsamples, "Samples &nbsp;&nbsp;&nbsp;", ncontrasts, "Comparisons")
        )
      } else {
        return(
          paste(ndatasets, "Data sets &nbsp;&nbsp;&nbsp;", nsamples, "Samples")
        )
      }
    })

    output$pgx_stats_ui <- shiny::renderUI(HTML(pgx_stats()))

    ## ================================================================================
    ## Data sets table
    ## ================================================================================

    ## reactive value for updating table
    touchtable <- shiny::reactiveVal(0)


    ## ------------------------------------------------
    ## Board return object
    ## ------------------------------------------------
    res <- list(
      is_data_loaded = is_data_loaded
    )
    return(res)
  })
}
