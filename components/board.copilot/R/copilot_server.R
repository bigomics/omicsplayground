#' Copilot Board Server
#'
#' Orchestrates the AI copilot: agent lifecycle, async chat streaming,
#' dataset/session management, and evidence panel wiring.
#'
#' @param id Module namespace id.
#' @param pgx Global PGX reactiveValues shared across all boards.
#' @param pgx_dir Path to the directory containing .pgx dataset files.
#' @param chat_dir Path to the directory for persisting chat sessions.
#' @param docs_dir Path to the directory for uploaded documents.
#' @param maxturns Maximum number of user turns per session.
#' @param tiers Character vector of copilot tier identifiers.
#' @param is_data_loaded Optional reactiveVal signalling dataset loads.

CopilotBoardServer <- function(id, pgx = NULL, pgx_dir = NULL,
                               chat_dir, docs_dir,
                               maxturns = Inf,
                               tiers = "copilot-default",
                               is_data_loaded = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    info("[CopilotBoard] server initialised —",
         " tiers=", paste(tiers, collapse = ","),
         " chat_dir=", chat_dir,
         " docs_dir=", docs_dir)

    ## --- Agent state ---
    copilot <- shiny::reactiveVal(NULL)
    n_turns <- shiny::reactiveVal(0)
    active_dataset_name <- shiny::reactiveVal(NULL)
    busy <- shiny::reactiveVal(FALSE)
    current_tier <- shiny::reactiveVal(tiers[1])
    session_dirty <- shiny::reactiveVal(FALSE)
    session_generation <- shiny::reactiveVal(0L)
    restoring <- shiny::reactiveVal(FALSE)
    restore_started_at <- shiny::reactiveVal(NULL)
    restore_async_invoked_at <- shiny::reactiveVal(NULL)
    restore_result <- shiny::reactiveVal(NULL)
    replay_started_at <- shiny::reactiveVal(NULL)
    replay_message_count <- shiny::reactiveVal(0L)
    pending_plot_requests <- shiny::reactiveVal(list())
    plot_build_inflight <- shiny::reactiveVal(FALSE)
    trace_enabled <- function() isTRUE(getOption("copilot.trace", FALSE))
    copilot_dbg <- function(...) {
      if (trace_enabled()) dbg(...)
      invisible(NULL)
    }
    info("[CopilotBoard] save flow=sync write + new_chat=recreate")

    restore_task <- ExtendedTask$new(function(session_id, session_dir, pgx_to_inject) {
      future_promise({
        restore_started <- Sys.time()
        worker_started_at <- as.numeric(restore_started)
        restored_agent <- omicsagentovi::ovi_restore(
          session_id  = session_id,
          session_dir = session_dir,
          pgx         = pgx_to_inject,
          restore_pgx = "never"
        )
        if (!is.null(pgx_to_inject)) {
          ## Keep dataset label + PGX object consistent even when manifest dataset differs.
          try(omicsagentovi::ovi_set_pgx(restored_agent, pgx_to_inject), silent = TRUE)
        }
        worker_finished <- Sys.time()
        restore_ms <- round(as.numeric(difftime(worker_finished, restore_started, units = "secs")) * 1000, 1)
        list(
          agent = restored_agent,
          restore_ms = restore_ms,
          restore_mode = "async",
          bound_current_pgx = !is.null(pgx_to_inject),
          worker_started_at = worker_started_at,
          worker_finished_at = as.numeric(worker_finished)
        )
      }, seed = TRUE)
    })

    ## --- Persistent chat store (lazy: no dir is created per-session
    ## until the first successful turn triggers session_save) ---
    chat_store <- omicsagentovi::SessionStore(session_dir = chat_dir)

    ## --- Flush the active session on browser close ---
    session$onSessionEnded(function() {
      wrapper <- shiny::isolate(copilot())
      if (!is.null(wrapper) && shiny::isolate(n_turns()) > 0) {
        payload <- tryCatch(
          omicsagentovi::session_collect_payload(chat_store, wrapper$agent),
          error = function(e) NULL
        )
        if (!is.null(payload)) {
          payload$manifest$turn_count <- shiny::isolate(n_turns())
          payload$manifest$copilot_save_generation <- as.integer(shiny::isolate(session_generation()))
          try(omicsagentovi::session_write_payload(chat_store, payload), silent = TRUE)
        }
      }
    })

    ## --- Populate tier selector once on first flush ---
    session$onFlushed(function() {
      tier_choices <- setNames(tiers, vapply(tiers, copilot_tier_label, character(1)))
      shiny::updateSelectInput(session, "tier",
                               choices  = tier_choices,
                               selected = tiers[1])
    }, once = TRUE)

    ## --- Reactive wrapper around global pgx for sub-modules ---
    ## Depend on pgx$name only — avoids cascade invalidation from unrelated PGX slot changes.
    local_pgx <- shiny::reactive({
      name <- pgx$name
      shiny::req(name)
      shiny::isolate({
        snapshot <- shiny::reactiveValuesToList(pgx)
        copilot_as_pgx(snapshot)
      })
    })

    ## --- Sub-module servers ---
    history_refresh <- shiny::reactiveVal(0)

    selected_dataset <- copilot_panel_datasets_server("datasets", pgx_dir = pgx_dir)
    history <- copilot_panel_history_server(
      "history",
      chat_store = chat_store,
      refresh_trigger = history_refresh
    )
    copilot_panel_docs_server("docs", docs_dir = docs_dir)

    evidence <- copilot_panel_evidence_server("evidence", local_pgx = local_pgx)

    build_plot_label <- function(plot_type) {
      label <- gsub("_", " ", plot_type)
      paste0(toupper(substring(label, 1, 1)), substring(label, 2))
    }

    clear_pending_plots <- function() {
      if (length(shiny::isolate(pending_plot_requests())) > 0L || isTRUE(shiny::isolate(plot_build_inflight()))) {
        copilot_dbg("[CopilotBoard] clear_pending_plots",
            " pending=", length(shiny::isolate(pending_plot_requests())),
            " inflight=", isTRUE(shiny::isolate(plot_build_inflight())))
      }
      pending_plot_requests(list())
      plot_build_inflight(FALSE)
      invisible(NULL)
    }

    plot_request_key <- function(plot_type, args) {
      relevant_args <- list(
        plot_type = plot_type %||% "",
        contrast = args$contrast %||% NULL,
        features = args$features %||% NULL,
        collection = args$collection %||% NULL
      )
      jsonlite::toJSON(relevant_args, auto_unbox = TRUE, null = "null")
    }

    persist_payload_sync <- function(payload, origin = "session_save") {
      turn_count <- payload$manifest$turn_count
      if (is.null(turn_count)) turn_count <- 0L
      save_generation <- payload$manifest$copilot_save_generation
      if (is.null(save_generation)) save_generation <- NA_integer_
      write_started_at <- Sys.time()
      omicsagentovi::session_write_payload(chat_store, payload)
      write_finished_at <- Sys.time()
      write_time_s <- as.numeric(difftime(write_finished_at, write_started_at, units = "secs"))
      list(
        session_id = payload$session_id,
        turn_count = turn_count,
        save_generation = as.integer(save_generation),
        write_time_s = write_time_s,
        origin = origin
      )
    }

    handle_saved_result <- function(result) {
      current_generation <- shiny::isolate(session_generation())
      dirty_before <- isTRUE(shiny::isolate(session_dirty()))
      copilot_dbg(
        "[CopilotBoard] saved session",
        " origin=", result$origin %||% "unknown",
        " session_id=", result$session_id %||% NA_character_,
        " turns=", result$turn_count %||% 0L,
        " save_generation=", result$save_generation %||% NA_integer_,
        " write_time_s=", round(as.numeric(result$write_time_s %||% NA_real_), 3),
        " current_generation=", current_generation,
        " dirty_before=", dirty_before
      )
      saved_generation <- result$save_generation %||% NA_integer_
      if (!is.na(saved_generation) && saved_generation >= current_generation) {
        session_dirty(FALSE)
      }
      copilot_dbg("[CopilotBoard] saved session dirty_after=", isTRUE(shiny::isolate(session_dirty())))
      tryCatch(
        copilot_prune_sessions(chat_store),
        error = function(e) info("[CopilotBoard] prune failed: ", conditionMessage(e))
      )
      history_refresh(shiny::isolate(history_refresh()) + 1L)
      invisible(NULL)
    }

    schedule_session_persist <- function(agent_wrapper, turns_used) {
      session$onFlushed(function() {
        if (!isTRUE(shiny::isolate(session_dirty()))) {
          copilot_dbg("[CopilotBoard] session_save onFlushed skipped dirty=FALSE")
          return(invisible(NULL))
        }
        session_id <- tryCatch(agent_wrapper$agent@context@session_id, error = function(e) NA_character_)
        copilot_dbg(
          "[CopilotBoard] session_save onFlushed start",
          " session_id=", session_id,
          " turns_used=", turns_used,
          " dirty=", isTRUE(shiny::isolate(session_dirty())),
          " generation=", shiny::isolate(session_generation())
        )
        persist_start <- Sys.time()
        result <- tryCatch({
          payload <- omicsagentovi::session_collect_payload(chat_store, agent_wrapper$agent)
          payload$manifest$turn_count <- turns_used
          payload$manifest$copilot_save_generation <- as.integer(shiny::isolate(session_generation()))
          persist_payload_sync(payload, origin = "session_save")
        }, error = function(e) {
          info("[CopilotBoard] session_save failed: ", conditionMessage(e))
          NULL
        })
        if (!is.null(result)) handle_saved_result(result)
        persist_end <- Sys.time()
        copilot_dbg(
          "[CopilotBoard] session_save",
          " total_time=", round(as.numeric(difftime(persist_end, persist_start, units = "secs")), 3), "s",
          " write_time=", round(as.numeric(result$write_time_s %||% NA_real_), 3), "s"
        )
      }, once = TRUE)
      invisible(NULL)
    }

    mark_session_dirty <- function() {
      next_generation <- shiny::isolate(session_generation()) + 1L
      session_generation(next_generation)
      session_dirty(TRUE)
      next_generation
    }

    queue_plot_request <- function(pgx, plot_type, args) {
      reqs <- shiny::isolate(pending_plot_requests())
      queued_at <- Sys.time()
      req_key <- plot_request_key(plot_type, args %||% list())
      if (length(reqs) > 0L && any(vapply(reqs, function(x) identical(x$key, req_key), logical(1L)))) {
        return(invisible(NULL))
      }
      reqs[[length(reqs) + 1L]] <- list(
        pgx       = copilot_as_pgx(pgx),
        plot_type = plot_type,
        args      = args,
        timestamp = queued_at,
        queued_at = queued_at,
        key       = req_key,
        label     = build_plot_label(plot_type)
      )
      pending_plot_requests(reqs)
      copilot_dbg(
        "[CopilotBoard] queued plot request",
        " plot_type=", plot_type,
        " pending=", length(reqs),
        " args=", paste(names(args %||% list()), collapse = ",")
      )
      invisible(NULL)
    }

    flush_pending_plots <- function() {
      if (isTRUE(shiny::isolate(plot_build_inflight()))) return(invisible(NULL))
      reqs <- shiny::isolate(pending_plot_requests())
      if (!length(reqs)) return(invisible(NULL))

      copilot_dbg("[CopilotBoard] about to flush queued plots pending=", length(reqs))
      pending_plot_requests(list())
      plot_build_inflight(TRUE)
      on.exit(plot_build_inflight(FALSE), add = TRUE)

      for (req in reqs) {
        build_started_at <- Sys.time()
        tryCatch({
          copilot_dbg("[CopilotBoard] build_plot start",
              " plot_type=", req$plot_type)
          plot_obj <- copilot_build_plot(req$pgx, req$plot_type, req$args)
          kind <- copilot_detect_plot_kind(plot_obj)
          if (is.null(kind)) {
            stop("Unsupported plot object for evidence panel.", call. = FALSE)
          }

          ## Pre-render outside flushReact: ggplot -> PNG, plotly -> build
          prerendered_path <- NULL
          if (identical(kind, "ggplot")) {
            prerendered_path <- copilot_prerender_ggplot(plot_obj)
          }
          if (identical(kind, "plotly")) {
            plot_obj <- copilot_prerender_plotly(plot_obj)
          }

          built_at <- Sys.time()
          copilot_dbg("[CopilotBoard] build_plot done",
              " plot_type=", req$plot_type,
              " build_time=",
              round(as.numeric(difftime(built_at, build_started_at, units = "secs")), 3), "s")

          evidence$append_plot(list(
            plot      = plot_obj,
            kind      = kind,
            prerendered_path = prerendered_path,
            plot_type = req$plot_type,
            args      = req$args,
            timestamp = req$timestamp,
            queued_at = req$queued_at,
            built_at  = built_at,
            label     = req$label
          ))
          copilot_dbg("[CopilotBoard] append_plot done plot_type=", req$plot_type)
        }, error = function(e) {
          info("[CopilotBoard] deferred plot build failed: ", conditionMessage(e))
          shiny::showNotification(
            paste0("Failed to render plot: ", conditionMessage(e)),
            type = "error",
            duration = 6
          )
        })
      }
      invisible(NULL)
    }

    ## --- Plot callback for evidence panel ---
    ## Queue plot requests during streaming and build them once the stream settles.
    plot_callback_impl <- function(pgx, plot_type, args) {
      queue_plot_request(pgx, plot_type, args)
    }
    plot_callback <- if (is.function(session$wrapFunction)) {
      session$wrapFunction(plot_callback_impl)
    } else {
      function(...) shiny::withReactiveDomain(session, plot_callback_impl(...))
    }

    ## --- Greeting (shown once before any dataset) ---
    session$onFlushed(function() {
      mesg <- "Hi, I'm **BigOmics Copilot**! Load a dataset to get started."
      shinychat::chat_append("chat", mesg)
    }, once = TRUE)

    ## --- Create/reset agent whenever the GLOBAL PGX changes ---
    ## This fires both when the user loads data via DataView/Home and when
    ## the copilot's own left panel writes into `pgx` below.
    shiny::observeEvent(pgx$name, {
      shiny::req(pgx$name)
      if (isTRUE(restoring())) return()
      dataset_name <- as.character(pgx$name[[1]])
      pgx_x <- shiny::isolate(pgx$X)
      shiny::req(pgx_x)
      pgx_snapshot <- copilot_snapshot_pgx(pgx)

      params <- list(
        name  = dataset_name,
        nsamp = omicsai::omicsai_format_num(ncol(pgx_x), digits = 0),
        ngene = omicsai::omicsai_format_num(nrow(pgx_x), digits = 0)
      )

      if (is.null(copilot())) {
        info("[CopilotBoard] first dataset — creating agent for: ", dataset_name)
        agent_wrapper <- copilot_create_agent(
          pgx = pgx_snapshot,
          plot_callback = plot_callback,
          pgx_dir = pgx_dir,
          docs_dir = docs_dir,
          tier = current_tier()
        )
        model_name <- tryCatch(agent_wrapper$agent@model, error = function(e) "unknown")
        info("[CopilotBoard] agent created — model=", model_name)
        copilot(agent_wrapper)
        active_dataset_name(dataset_name)
        n_turns(0)
        session_generation(0L)
        session_dirty(FALSE)
        evidence$clear_plots()
        evidence$clear_table()
        shinychat::chat_clear("chat")
        shinychat::chat_append("chat", omicsai::omicsai_substitute_template(
          "Dataset **{{name}}** loaded — {{nsamp}} samples, {{ngene}} genes. What would you like to know?",
          params
        ))
      } else if (!identical(active_dataset_name(), dataset_name)) {
        info("[CopilotBoard] dataset switch in place: ",
             active_dataset_name(), " -> ", dataset_name)
        copilot()$set_pgx(pgx_snapshot)
        active_dataset_name(dataset_name)
        session_dirty(FALSE)
        shinychat::chat_append("chat", omicsai::omicsai_substitute_template(
          "_Switched to dataset **{{name}}** — {{nsamp}} samples, {{ngene}} genes._",
          params
        ))
      }
    }, ignoreNULL = TRUE)

    ## --- When user picks a dataset in the copilot's left panel, load it
    ## into the GLOBAL PGX so it propagates to all boards. ---
    shiny::observeEvent(selected_dataset(), {
      dataset_path <- selected_dataset()
      shiny::req(dataset_path)
      shiny::req(file.exists(dataset_path))

      ## Avoid reloading if already the active dataset
      dataset_name <- sub("[.]pgx$", "", basename(dataset_path))
      if (!is.null(pgx$name) && pgx$name == dataset_name) {
        info("[CopilotBoard] dataset already active, skipping reload: ", dataset_name)
        return()
      }

      info("[CopilotBoard] loading PGX into global state: ", dataset_path)
      load_result <- shiny::withProgress(message = "Loading dataset for Copilot...", {
        tryCatch({
          p <- playbase::pgx.load(dataset_path)
          p <- playbase::pgx.initialize(p)
          ## Post-condition: pgx.initialize silently returns NULL for
          ## non-PGX lists, so validate shape before trusting it.
          if (is.null(p) || !is.list(p) || is.null(p$X)) {
            stop("file did not yield a valid PGX object (missing expression matrix)")
          }
          p$name <- dataset_name
          list(ok = TRUE, pgx = p)
        }, error = function(e) {
          list(ok = FALSE, err = conditionMessage(e))
        })
      })

      if (!load_result$ok) {
        warning("[CopilotBoard] failed to load PGX: ", dataset_path,
                " :: ", load_result$err)
        shiny::showNotification(
          paste0("Failed to load dataset '", dataset_name, "': ", load_result$err),
          type = "error",
          duration = 8
        )
        return()
      }
      loaded_pgx <- load_result$pgx

      ## Copy into global pgx reactiveValues (same semantics as loading_server.R)
      shiny::isolate(sync_rv_from_list(pgx, loaded_pgx))
      info("[CopilotBoard] PGX loaded into global — samples=", ncol(pgx$X), " genes=", nrow(pgx$X))

      ## Signal the app that a new dataset was loaded so sidebar/tabs refresh
      if (is.function(is_data_loaded)) {
        val <- is_data_loaded()
        is_data_loaded(if (is.null(val)) 1L else val + 1L)
      }

      ## Refresh sidebar menu visibility. The is_data_loaded observer in server.R
      ## also calls this, but the direct call here ensures the JS fires immediately
      ## without waiting for the reactive flush cycle.
      bigdash.showTabsGoToDataView(session)
    })

    ## --- Core ask function ---
    ask_copilot <- function(question, showq = TRUE) {
      agent_wrapper <- copilot()
      if (is.null(agent_wrapper)) {
        shinychat::chat_append("chat", "Please load a dataset first.")
        return(NULL)
      }
      if (n_turns() >= maxturns) {
        shinychat::chat_append("chat",
          paste("You've reached the", maxturns, "turn limit. Click **+ New chat** to start a new conversation."))
        return(NULL)
      }
      if (busy()) {
        shinychat::chat_append("chat", "Please wait for the current response to finish.")
        return(NULL)
      }
      if (showq) {
        msg <- list(role = "user", content = question)
        shinychat::chat_append_message("chat", msg, chunk = FALSE)
      }

      busy(TRUE)
      copilot_dbg("[CopilotBoard] ask_copilot start question_nchar=", nchar(question))

      ## Async path: stream_async returns a coro generator of Content objects.
      ## shinychat::chat_append drives the generator via later(), yielding
      ## each ContentText / ContentToolRequest / ContentToolResult to the
      ## browser as it arrives — no manual callbacks needed.
      tryCatch({
        stream <- agent_wrapper$stream_async(question)
        result <- shinychat::chat_append("chat", stream)

        promises::then(result,
          onFulfilled = function(value) {
            copilot_dbg("[CopilotBoard] stream fulfilled")
            turns_used <- shiny::isolate(n_turns()) + 1L
            n_turns(turns_used)
            mark_session_dirty()
            tryCatch({
              had_plots <- length(shiny::isolate(pending_plot_requests())) > 0L
              flush_pending_plots()
              if (had_plots && is.function(evidence$plot_rendered)) {
                copilot_dbg("[CopilotBoard] plot_rendered — scheduling session_save")
                shiny::observeEvent(evidence$plot_rendered(), {
                  schedule_session_persist(agent_wrapper, turns_used = turns_used)
                }, once = TRUE, ignoreInit = TRUE)
              } else {
                later::later(function() {
                  schedule_session_persist(agent_wrapper, turns_used = turns_used)
                }, delay = 0.5)
              }
            }, error = function(e) {
              info("[CopilotBoard] session_save failed: ", conditionMessage(e))
            })
          },
          onRejected = function(err) {
            clear_pending_plots()
            info("[CopilotBoard] prompt failed: ", conditionMessage(err))
            shinychat::chat_append("chat", paste("Copilot error:", conditionMessage(err)))
          }
        ) |> promises::finally(function() {
          copilot_dbg("[CopilotBoard] ask_copilot finally busy=FALSE")
          busy(FALSE)
        })
      }, error = function(e) {
        clear_pending_plots()
        info("[CopilotBoard] stream_async failed: ", conditionMessage(e))
        shinychat::chat_append("chat", paste("Copilot error:", conditionMessage(e)))
        busy(FALSE)
      })
    }

    ## --- User chat input ---
    shiny::observeEvent(input$chat_user_input, {
      shiny::req(input$chat_user_input)
      ask_copilot(input$chat_user_input, showq = FALSE)
    })

    ## --- Tier change ---
    ## Changing tier rebuilds the agent on the new (model, reasoning_effort)
    ## and resets the chat. History is not preserved in v1 — keeps the code
    ## simple and matches the existing Reset UX.
    shiny::observeEvent(input$tier, {
      new_tier <- input$tier
      if (is.null(new_tier) || !nzchar(new_tier)) return()
      if (identical(new_tier, current_tier())) return()
      if (busy()) {
        shinychat::chat_append("chat",
          "Please wait for the current response to finish before switching tiers.")
        shiny::updateSelectInput(session, "tier", selected = current_tier())
        return()
      }
      info("[CopilotBoard] tier change — ", current_tier(), " -> ", new_tier)
      current_tier(new_tier)

      if (!is.null(pgx$X)) {
        agent_wrapper <- copilot_create_agent(
          pgx = copilot_snapshot_pgx(pgx),
          plot_callback = plot_callback,
          pgx_dir = pgx_dir,
          docs_dir = docs_dir,
          tier = new_tier
        )
        copilot(agent_wrapper)
      }
      n_turns(0)
      session_generation(0L)
      session_dirty(FALSE)
      clear_pending_plots()
      evidence$clear_plots()
      evidence$clear_table()
      shinychat::chat_clear("chat")
      shinychat::chat_append("chat",
        paste0("Switched to **", copilot_tier_label(new_tier), "** tier. Chat reset."))
      shiny::showNotification(
        paste0("Copilot tier: ", copilot_tier_label(new_tier)),
        type = "message", duration = 3
      )
    }, ignoreInit = TRUE)

    ## --- New chat ---
    shiny::observeEvent(input$new_chat, {
      new_chat_started_at <- Sys.time()
      if (is.null(input$new_chat) || input$new_chat < 1) return()
      if (busy()) {
        shinychat::chat_append("chat", "Please wait for the current response to finish.")
        return()
      }
      ## Flush current session synchronously before resetting.
      wrapper <- copilot()
      needs_preflush <- !is.null(wrapper) && n_turns() > 0 && isTRUE(session_dirty())
      if (needs_preflush) {
        preflush_session_id <- tryCatch(wrapper$agent@context@session_id, error = function(e) NA_character_)
        copilot_dbg(
          "[CopilotBoard] new chat preflush start",
          " session_id=", preflush_session_id,
          " n_turns=", n_turns(),
          " dirty=", isTRUE(session_dirty())
        )
        preflush_started_at <- Sys.time()
        result <- tryCatch({
          payload <- omicsagentovi::session_collect_payload(chat_store, wrapper$agent)
          persist_payload_sync(payload, origin = "new_chat_preflush")
        }, error = function(e) {
          info("[CopilotBoard] new chat preflush save failed: ", conditionMessage(e))
          NULL
        })
        if (!is.null(result)) handle_saved_result(result)
        preflush_elapsed <- round(as.numeric(difftime(Sys.time(), preflush_started_at, units = "secs")) * 1000, 1)
        copilot_dbg("[CopilotBoard] new chat preflush done elapsed_ms=", preflush_elapsed)
      } else {
        copilot_dbg("[CopilotBoard] new chat preflush skipped dirty=", isTRUE(session_dirty()))
      }
      info("[CopilotBoard] new chat requested — turns so far=", n_turns())
      if (!is.null(pgx$X)) {
        wrapper <- copilot()
        old_session_id <- tryCatch(wrapper$agent@context@session_id, error = function(e) NA_character_)
        snapshot_started_at <- Sys.time()
        pgx_snapshot <- copilot_snapshot_pgx(pgx)
        snapshot_elapsed_ms <- round(as.numeric(difftime(Sys.time(), snapshot_started_at, units = "secs")) * 1000, 1)

        create_started_at <- Sys.time()
        agent_wrapper <- copilot_create_agent(
          pgx = pgx_snapshot,
          plot_callback = plot_callback,
          pgx_dir = pgx_dir,
          docs_dir = docs_dir,
          tier = current_tier()
        )
        create_elapsed_ms <- round(as.numeric(difftime(Sys.time(), create_started_at, units = "secs")) * 1000, 1)
        new_session_id <- tryCatch(agent_wrapper$agent@context@session_id, error = function(e) NA_character_)
        info(
          "[CopilotBoard] new chat recreate",
          " snapshot_ms=", snapshot_elapsed_ms,
          " create_ms=", create_elapsed_ms,
          " old_session_id=", old_session_id,
          " new_session_id=", new_session_id
        )
        copilot(agent_wrapper)
        active_dataset_name(as.character(pgx$name[[1]]))
      } else {
        copilot(NULL)
        active_dataset_name(NULL)
      }
      n_turns(0)
      session_generation(0L)
      session_dirty(FALSE)
      clear_pending_plots()
      evidence$clear_plots()
      evidence$clear_table()
      shinychat::chat_clear("chat")
      shinychat::chat_append("chat", "New conversation started. Ask me anything about your data!")
      info(
        "[CopilotBoard] new chat done total_ms=",
        round(as.numeric(difftime(Sys.time(), new_chat_started_at, units = "secs")) * 1000, 1)
      )
    }, ignoreInit = TRUE)

    ## --- Restore session (async) ---
    emit_restore_result <- function(payload) {
      restore_result(list(payload = payload, ts = as.numeric(Sys.time())))
      invisible(NULL)
    }

    complete_restore <- function(restored_payload) {
      if (is.null(restored_payload)) {
        shinychat::chat_clear("chat")
        shinychat::chat_append("chat", "Failed to restore session.")
        restore_started_at(NULL)
        restoring(FALSE)
        return(invisible(NULL))
      }
      restored <- restored_payload$agent %||% restored_payload
      restore_ms <- restored_payload$restore_ms %||% NA_real_
      restore_mode <- restored_payload$restore_mode %||% "unknown"
      async_roundtrip_ms <- NA_real_
      async_queue_ms <- NA_real_
      async_worker_wall_ms <- NA_real_
      async_post_ms <- NA_real_
      if (identical(restore_mode, "async")) {
        invoked_at <- tryCatch(shiny::isolate(restore_async_invoked_at()), error = function(e) NULL)
        completed_at <- Sys.time()
        if (!is.null(invoked_at)) {
          async_roundtrip_ms <- round(as.numeric(difftime(completed_at, invoked_at, units = "secs")) * 1000, 1)
        }
        worker_started_at <- restored_payload$worker_started_at %||% NA_real_
        worker_finished_at <- restored_payload$worker_finished_at %||% NA_real_
        if (!is.na(worker_started_at) && !is.na(worker_finished_at)) {
          async_worker_wall_ms <- round((worker_finished_at - worker_started_at) * 1000, 1)
        }
        if (!is.null(invoked_at) && !is.na(worker_started_at)) {
          async_queue_ms <- round((worker_started_at - as.numeric(invoked_at)) * 1000, 1)
        }
        if (!is.na(async_roundtrip_ms) && !is.na(async_queue_ms) && !is.na(async_worker_wall_ms)) {
          async_post_ms <- round(async_roundtrip_ms - async_queue_ms - async_worker_wall_ms, 1)
          if (async_post_ms < 0) async_post_ms <- 0
        }
      }
      restore_async_invoked_at(NULL)

      ## Bind copilot-side state onto the restored agent
      restored_agent <- restored
      bound_current_pgx <- isTRUE(restored_payload$bound_current_pgx %||% FALSE)
      if (is.function(plot_callback)) {
        restored_agent@context@state$plot_callback <- plot_callback
      }
      if (!is.null(pgx_dir)) {
        restored_agent@context@state$data_dir <- pgx_dir
      }
      if (!is.null(docs_dir)) {
        restored_agent@context@state$docs_dir <- docs_dir
      }
      copilot_dbg(
        "[CopilotBoard] restored agent rebound",
        " has_plot_callback=", is.function(restored_agent@context@state$plot_callback),
        " has_data_dir=", !is.null(restored_agent@context@state$data_dir)
      )
      current_tier_value <- tryCatch(shiny::isolate(current_tier()), error = function(e) tiers[1])
      wrapper <- list(
        prompt = function(text) omicsagentovi::agent_prompt(restored_agent, text),
        prompt_stream = function(text, on_event = NULL) {
          omicsagentovi::agent_prompt_stream(restored_agent, text, on_event = on_event)
        },
        stream_async = function(text) {
          stream <- omicsagentovi::agent_prompt_async(restored_agent, text)
          coro::async_generator(function() {
            for (item in coro::await_each(stream)) {
              if (S7::S7_inherits(item, ellmer::ContentToolRequest)) {
                coro::yield(ellmer::ContentText(.format_tool_request(item)))
              } else if (!S7::S7_inherits(item, ellmer::ContentToolResult)) {
                coro::yield(item)
              }
            }
          })()
        },
        agent = restored_agent,
        tier = current_tier_value,
        set_pgx = function(new_pgx) {
          omicsagentovi::ovi_set_pgx(restored_agent, copilot_as_pgx(new_pgx))
        },
        reset = function(...) invisible(NULL)
      )

      copilot(wrapper)
      clear_pending_plots()
      evidence$clear_plots()
      evidence$clear_table()
      shinychat::chat_clear("chat")

      turns <- tryCatch(restored_agent@chat$get_turns(), error = function(e) list())
      replay_collect_started_at <- Sys.time()
      replay_messages <- copilot_collect_replay_messages(turns, policy = "expected")
      replay_collect_ms <- round(as.numeric(difftime(Sys.time(), replay_collect_started_at, units = "secs")) * 1000, 1)
      replay_msg_count <- length(replay_messages)
      replay_message_count(replay_msg_count)
      replay_msg_nchars <- vapply(replay_messages, function(msg) nchar(msg$content %||% "", type = "bytes"), integer(1L))
      replay_chars_total <- sum(replay_msg_nchars)
      replay_chars_max <- if (length(replay_msg_nchars) > 0L) max(replay_msg_nchars) else 0L
      replay_mode <- getOption("copilot.restore_replay_mode", "single")
      if (!identical(replay_mode, "batch")) replay_mode <- "single"

      replay_apply_started_at <- Sys.time()
      n_user <- copilot_replay_turns(
        "chat",
        turns,
        policy = "expected",
        replay_messages = replay_messages,
        mode = replay_mode,
        session = session,
        batch_size = 32L,
        done_input_id = session$ns("chat_replay_done")
      )
      replay_apply_ms <- round(as.numeric(difftime(Sys.time(), replay_apply_started_at, units = "secs")) * 1000, 1)
      n_turns(n_user)
      session_generation(0L)
      session_dirty(FALSE)
      active_dataset_name(NULL)
      history$restore_request(NULL)
      restoring(FALSE)
      total_restore_ms <- NA_real_
      restore_started_snapshot <- tryCatch(shiny::isolate(restore_started_at()), error = function(e) NULL)
      if (!is.null(restore_started_snapshot)) {
        total_restore_ms <- round(as.numeric(difftime(Sys.time(), restore_started_snapshot, units = "secs")) * 1000, 1)
      }
      restore_started_at(NULL)
      info(
        "[CopilotBoard] session restored — user turns replayed=", n_user,
        " restore_mode=", restore_mode,
        " restore_ms=", restore_ms,
        " total_ms=", total_restore_ms,
        " async_roundtrip_ms=", async_roundtrip_ms,
        " async_queue_ms=", async_queue_ms,
        " async_worker_wall_ms=", async_worker_wall_ms,
        " async_post_ms=", async_post_ms,
        " bound_current_pgx=", bound_current_pgx,
        " replay_mode=", replay_mode,
        " replay_messages=", replay_msg_count,
        " replay_chars_total=", replay_chars_total,
        " replay_chars_max=", replay_chars_max,
        " replay_collect_ms=", replay_collect_ms,
        " replay_apply_ms=", replay_apply_ms
      )

      if (!identical(replay_mode, "batch")) {
        replay_started_at(NULL)
        replay_message_count(0L)
      } else {
        replay_started_at(replay_apply_started_at)
      }
      invisible(NULL)
    }

    shiny::observeEvent(history$restore_request(), {
      session_id <- history$restore_request()
      shiny::req(session_id)
      if (busy()) {
        shinychat::chat_append("chat", "Please wait for the current response to finish.")
        return()
      }

      ## Conversation-only restore policy: never reload PGX from disk.
      ## Always bind to the currently loaded app dataset when available.
      pgx_to_inject <- tryCatch(local_pgx(), error = function(e) NULL)
      if (!is.null(pgx_to_inject)) {
        copilot_dbg("[CopilotBoard] restore policy — binding current app PGX: ", as.character(pgx_to_inject$name %||% "unknown"))
      }

      info("[CopilotBoard] restoring session: ", session_id)
      restoring(TRUE)
      restore_started_at(Sys.time())
      shinychat::chat_clear("chat")
      shinychat::chat_append("chat", "Restoring session...")

      if (!is.null(pgx_to_inject)) {
        restore_async_invoked_at(NULL)
        session$onFlushed(function() {
          sync_started <- Sys.time()
          restored_payload <- tryCatch({
            restored_agent <- omicsagentovi::ovi_restore(
              session_id  = session_id,
              session_dir = chat_store@state$session_dir,
              pgx         = pgx_to_inject,
              restore_pgx = "never"
            )
            if (!is.null(pgx_to_inject)) {
              ## Keep dataset label + PGX object consistent even when manifest dataset differs.
              try(omicsagentovi::ovi_set_pgx(restored_agent, pgx_to_inject), silent = TRUE)
            }
            restore_ms <- round(as.numeric(difftime(Sys.time(), sync_started, units = "secs")) * 1000, 1)
            list(
              agent = restored_agent,
              restore_ms = restore_ms,
              restore_mode = "sync",
              bound_current_pgx = !is.null(pgx_to_inject)
            )
          }, error = function(e) {
            info("[CopilotBoard] restore sync failed: ", conditionMessage(e))
            NULL
          })
          emit_restore_result(restored_payload)
        }, once = TRUE)
        return()
      }

      tryCatch(
        {
          restore_async_invoked_at(Sys.time())
          restore_task$invoke(session_id, chat_store@state$session_dir, pgx_to_inject)
        },
        error = function(e) {
          info("[CopilotBoard] restore_task invoke failed: ", conditionMessage(e))
          shinychat::chat_clear("chat")
          shinychat::chat_append("chat", "Failed to restore session.")
          restore_started_at(NULL)
          restore_async_invoked_at(NULL)
          restoring(FALSE)
        }
      )
    }, ignoreNULL = TRUE)

    ## --- Handle restore completion ---
    shiny::observeEvent(restore_task$result(), {
      restored_payload <- tryCatch(restore_task$result(), error = function(e) NULL)
      emit_restore_result(restored_payload)
    }, ignoreNULL = TRUE)

    shiny::observeEvent(restore_result(), {
      event <- restore_result()
      shiny::req(event)
      restore_result(NULL)
      complete_restore(event$payload)
    }, ignoreNULL = TRUE)

    shiny::observeEvent(input$chat_replay_done, {
      started_at <- replay_started_at()
      if (is.null(started_at)) return()
      replay_ms <- round(as.numeric(difftime(Sys.time(), started_at, units = "secs")) * 1000, 1)
      info(
        "[CopilotBoard] replay complete mode=batch",
        " replay_messages=", replay_message_count(),
        " replay_ms=", replay_ms
      )
      replay_started_at(NULL)
      replay_message_count(0L)
    }, ignoreNULL = TRUE, ignoreInit = TRUE)

    ## --- Delete session ---
    shiny::observeEvent(history$delete_request(), {
      session_id <- history$delete_request()
      shiny::req(session_id)
      info("[CopilotBoard] deleting session: ", session_id)
      tryCatch(
        omicsagentovi::ovi_session_delete(session_id, session_dir = chat_store@state$session_dir),
        error = function(e) warning("[CopilotBoard] session delete failed: ", conditionMessage(e))
      )
      history$delete_request(NULL)
      history_refresh(history_refresh() + 1L)
    }, ignoreNULL = TRUE)

    ## --- Example question buttons ---
    shiny::observeEvent(input$ask_describe, {
      shiny::req(copilot())
      ask_copilot("Describe my experiment and the main comparisons being made")
    })
    shiny::observeEvent(input$ask_findings, {
      shiny::req(copilot())
      ask_copilot("Summarize the main findings of this experiment")
    })
    shiny::observeEvent(input$ask_pathways, {
      shiny::req(copilot())
      ask_copilot("List the top enriched pathways and explain how they relate to the experiment")
    })
    shiny::observeEvent(input$ask_biomarkers, {
      shiny::req(copilot())
      ask_copilot("Show the top candidate biomarkers for this experiment")
    })
    shiny::observeEvent(input$ask_plot, {
      shiny::req(copilot())
      ask_copilot("Show me a PCA plot of the samples")
    })
  })
}
