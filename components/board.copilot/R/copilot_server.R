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
                               maxturns = 10,
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
    restoring <- shiny::reactiveVal(FALSE)
    pending_plot_requests <- shiny::reactiveVal(list())
    plot_build_inflight <- shiny::reactiveVal(FALSE)

    ## --- Session save backend ---
    save_task <- ExtendedTask$new(function(session_dir, payload) {
      future_promise({
        turn_count <- payload$manifest$turn_count
        if (is.null(turn_count)) turn_count <- 0L
        store <- omicsagentovi::SessionStore(session_dir = session_dir)
        omicsagentovi::session_write_payload(store, payload)
        list(session_id = payload$session_id, turn_count = turn_count)
      })
    })

    restore_task <- ExtendedTask$new(function(session_id, session_dir, pgx_to_inject) {
      future_promise({
        omicsagentovi::ovi_restore(
          session_id  = session_id,
          session_dir = session_dir,
          pgx         = pgx_to_inject
        )
      })
    })

    ## --- Persistent chat store (lazy: no dir is created per-session
    ## until the first successful turn triggers session_save) ---
    chat_store <- omicsagentovi::SessionStore(session_dir = chat_dir)

    ## --- Flush the active session on browser close ---
    session$onSessionEnded(function() {
      wrapper <- shiny::isolate(copilot())
      if (!is.null(wrapper) && shiny::isolate(n_turns()) > 0) {
        ## Async save — avoid blocking session teardown with a 10-20s serialize.
        payload <- tryCatch(
          omicsagentovi::session_collect_payload(chat_store, wrapper$agent),
          error = function(e) NULL
        )
        if (!is.null(payload)) {
          tryCatch(
            save_task$invoke(chat_store@state$session_dir, payload),
            error = function(e) {
              ## save_task may already be running; fall back to sync.
              try(omicsagentovi::session_save(chat_store, wrapper$agent), silent = TRUE)
            }
          )
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
      if (length(pending_plot_requests()) > 0L || isTRUE(plot_build_inflight())) {
        if (getOption("copilot.trace", FALSE)) {
          dbg("[CopilotBoard] clear_pending_plots",
              " pending=", length(pending_plot_requests()),
              " inflight=", isTRUE(plot_build_inflight()))
        }
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

    schedule_session_persist <- function(agent_wrapper, turns_used) {
      ## Drop intermediate saves — only the latest session state matters.
      ## ExtendedTask queues invocations, so stale saves would pile up
      ## and each would run session_collect_payload synchronously in onFlushed.
      if (identical(shiny::isolate(save_task$status()), "running")) {
        dbg("[CopilotBoard] save already pending — skipping redundant persist")
        return(invisible(NULL))
      }
      session$onFlushed(function() {
        collect_start <- Sys.time()
        if (getOption("copilot.trace", FALSE)) {
          dbg("[CopilotBoard] session_save collect_start")
        }
        payload <- omicsagentovi::session_collect_payload(chat_store, agent_wrapper$agent)
        payload$manifest$turn_count <- turns_used
        collect_end <- Sys.time()

        invoke_start <- Sys.time()
        if (getOption("copilot.trace", FALSE)) {
          dbg("[CopilotBoard] session_save invoke_start")
        }
        save_task$invoke(chat_store@state$session_dir, payload)
        invoke_end <- Sys.time()

        dbg("[CopilotBoard] session_save",
            " collect_time=",
            round(as.numeric(difftime(collect_end, collect_start, units = "secs")), 3), "s",
            " invoke_time=",
            round(as.numeric(difftime(invoke_end, invoke_start, units = "secs")), 3), "s")
      }, once = TRUE)
      invisible(NULL)
    }

    queue_plot_request <- function(pgx, plot_type, args) {
      reqs <- pending_plot_requests()
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
      dbg(
        "[CopilotBoard] queued plot request",
        " plot_type=", plot_type,
        " pending=", length(reqs),
        " args=", paste(names(args %||% list()), collapse = ",")
      )
      invisible(NULL)
    }

    flush_pending_plots <- function() {
      if (isTRUE(plot_build_inflight())) return(invisible(NULL))
      reqs <- pending_plot_requests()
      if (!length(reqs)) return(invisible(NULL))

      if (getOption("copilot.trace", FALSE)) {
        dbg("[CopilotBoard] about to flush queued plots pending=", length(reqs))
      }
      pending_plot_requests(list())
      plot_build_inflight(TRUE)
      on.exit(plot_build_inflight(FALSE), add = TRUE)

      for (req in reqs) {
        build_started_at <- Sys.time()
        tryCatch({
          dbg("[CopilotBoard] build_plot start",
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
          dbg("[CopilotBoard] build_plot done",
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
          if (getOption("copilot.trace", FALSE)) {
            dbg("[CopilotBoard] append_plot done plot_type=", req$plot_type)
          }
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

    shiny::observeEvent(save_task$result(), {
      result <- save_task$result()
      shiny::req(result)
      info("[CopilotBoard] saved session turns=", result$turn_count %||% 0L)
      tryCatch(
        copilot_prune_sessions(chat_store),
        error = function(e) info("[CopilotBoard] prune failed: ", conditionMessage(e))
      )
      history_refresh(history_refresh() + 1L)
    }, ignoreNULL = TRUE)

    shiny::observeEvent(save_task$status(), {
      status <- save_task$status()
      if (!identical(status, "error")) return()
      err <- tryCatch(save_task$result(), error = function(e) e)
      info("[CopilotBoard] session_save failed: ", conditionMessage(err))
    }, ignoreInit = TRUE)

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
          tier = current_tier()
        )
        model_name <- tryCatch(agent_wrapper$agent@model, error = function(e) "unknown")
        info("[CopilotBoard] agent created — model=", model_name)
        copilot(agent_wrapper)
        active_dataset_name(dataset_name)
        n_turns(0)
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
        info("[CopilotBoard] ask_copilot called but no agent — no dataset loaded")
        shinychat::chat_append("chat", "Please load a dataset first.")
        return(NULL)
      }
      if (n_turns() >= maxturns) {
        info("[CopilotBoard] turn limit reached — maxturns=", maxturns)
        shinychat::chat_append("chat",
          paste("You've reached the", maxturns, "turn limit. Click **+ New chat** in the History panel to start a new conversation."))
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
      dbg("[CopilotBoard] ask_copilot start question_nchar=", nchar(question))

      ## Async path: stream_async returns a coro generator of Content objects.
      ## shinychat::chat_append drives the generator via later(), yielding
      ## each ContentText / ContentToolRequest / ContentToolResult to the
      ## browser as it arrives — no manual callbacks needed.
      tryCatch({
        stream <- agent_wrapper$stream_async(question)
        result <- shinychat::chat_append("chat", stream)

        promises::then(result,
          onFulfilled = function(value) {
            dbg("[CopilotBoard] stream fulfilled")
            turns_used <- n_turns() + 1L
            n_turns(turns_used)
            tryCatch({
              had_plots <- length(pending_plot_requests()) > 0L
              flush_pending_plots()
              if (had_plots && is.function(evidence$plot_rendered)) {
                dbg("[CopilotBoard] plot_rendered — scheduling session_save")
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
          dbg("[CopilotBoard] ask_copilot finally busy=FALSE")
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
          tier = new_tier
        )
        copilot(agent_wrapper)
      }
      n_turns(0)
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

    ## --- New chat (from History panel) ---
    shiny::observeEvent(history$new_chat_request(), {
      if (history$new_chat_request() == 0) return()
      if (busy()) {
        shinychat::chat_append("chat", "Please wait for the current response to finish.")
        return()
      }
      ## Flush current session asynchronously before resetting.
      wrapper <- copilot()
      if (!is.null(wrapper) && n_turns() > 0) {
        payload <- tryCatch(
          omicsagentovi::session_collect_payload(chat_store, wrapper$agent),
          error = function(e) NULL
        )
        if (!is.null(payload)) {
          tryCatch(
            save_task$invoke(chat_store@state$session_dir, payload),
            error = function(e) {
              info("[CopilotBoard] save_task busy on new-chat flush, dropping")
            }
          )
        }
      }
      info("[CopilotBoard] new chat requested — turns so far=", n_turns())
      if (!is.null(pgx$X)) {
        agent_wrapper <- copilot_create_agent(
          pgx = copilot_snapshot_pgx(pgx),
          plot_callback = plot_callback,
          pgx_dir = pgx_dir,
          tier = current_tier()
        )
        copilot(agent_wrapper)
        active_dataset_name(as.character(pgx$name[[1]]))
      } else {
        copilot(NULL)
        active_dataset_name(NULL)
      }
      n_turns(0)
      clear_pending_plots()
      evidence$clear_plots()
      evidence$clear_table()
      shinychat::chat_clear("chat")
      shinychat::chat_append("chat", "New conversation started. Ask me anything about your data!")
      history_refresh(history_refresh() + 1L)
    }, ignoreInit = TRUE)

    ## --- Restore session (async) ---
    shiny::observeEvent(history$restore_request(), {
      session_id <- history$restore_request()
      shiny::req(session_id)
      if (busy()) {
        shinychat::chat_append("chat", "Please wait for the current response to finish.")
        return()
      }

      ## Fast path: reuse already-loaded PGX when dataset matches
      pgx_to_inject <- NULL
      meta <- tryCatch(
        omicsagentovi::ovi_session_meta(session_id, session_dir = chat_store@state$session_dir),
        error = function(e) NULL
      )
      if (!is.null(meta) && !is.null(meta$dataset_name)) {
        current_pgx <- tryCatch(local_pgx(), error = function(e) NULL)
        if (!is.null(current_pgx) && identical(current_pgx$name, meta$dataset_name)) {
          pgx_to_inject <- current_pgx
          dbg("[CopilotBoard] restore fast path — reusing loaded PGX: ", meta$dataset_name)
        }
      }

      info("[CopilotBoard] restoring session: ", session_id)
      restoring(TRUE)
      shinychat::chat_clear("chat")
      shinychat::chat_append("chat", "Restoring session...")

      tryCatch(
        restore_task$invoke(session_id, chat_store@state$session_dir, pgx_to_inject),
        error = function(e) {
          info("[CopilotBoard] restore_task invoke failed: ", conditionMessage(e))
          shinychat::chat_clear("chat")
          shinychat::chat_append("chat", "Failed to restore session.")
          restoring(FALSE)
        }
      )
    }, ignoreNULL = TRUE)

    ## --- Handle restore completion ---
    shiny::observeEvent(restore_task$result(), {
      restored <- tryCatch(restore_task$result(), error = function(e) NULL)
      if (is.null(restored)) {
        shinychat::chat_clear("chat")
        shinychat::chat_append("chat", "Failed to restore session.")
        restoring(FALSE)
        return()
      }

      ## Bind copilot-side state onto the restored agent
      restored_agent <- restored
      if (is.function(plot_callback)) {
        restored_agent@context@state$plot_callback <- plot_callback
      }
      if (!is.null(pgx_dir)) {
        restored_agent@context@state$data_dir <- pgx_dir
      }
      dbg(
        "[CopilotBoard] restored agent rebound",
        " has_plot_callback=", is.function(restored_agent@context@state$plot_callback),
        " has_data_dir=", !is.null(restored_agent@context@state$data_dir)
      )
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
        tier = current_tier(),
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
      n_user <- copilot_replay_turns("chat", turns)
      n_turns(n_user)
      active_dataset_name(NULL)
      history$restore_request(NULL)
      restoring(FALSE)
      info("[CopilotBoard] session restored — user turns replayed=", n_user)
    }, ignoreNULL = TRUE)

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
