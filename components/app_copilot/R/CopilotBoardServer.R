# CopilotBoardServer.R — Copilot Board Orchestrator
#
# Lifecycle owner for the Agent reactive. Wires module servers and lifecycle
# controllers. Does not implement run/save/restore logic — those live in their
# respective controllers (copilot_run_controller, copilot_save_controller,
# copilot_restore_controller).

#' Copilot Board Orchestrator
#'
#' Wires module servers and controllers. Holds the agent reactive.
#'
#' @param id         Module namespace id.
#' @param pgx        Global PGX reactiveValues shared across all boards.
#' @param pgx_dir    Path to directory containing .pgx dataset files.
#' @param auth       Auth module reactiveValues (must expose `$user_dir` and
#'   `$email`). Copilot's `chat_dir` and `docs_dir` are always scoped to the
#'   per-user folder `<pgx_dir>/<email>` when an email is known, independent
#'   of ENABLE_USERDIR (see the path-resolution block below); with no email
#'   it falls back to `auth$user_dir`.
#' @param maxturns   Maximum user turns per session.
#' @param tiers      Character vector of tier identifiers (first = default).
#'   Defaults to `COPILOT_TIERS` from `copilot_options.R`.
#' @param is_data_loaded Optional reactiveVal; incremented after dataset load.
#'
#' @return NULL (side-effects only; called for Shiny registration).
#' @export
CopilotBoardServer <- function(
  id,
  pgx            = NULL,
  pgx_dir        = NULL,
  auth,
  maxturns       = Inf,
  tiers          = COPILOT_TIERS,
  is_data_loaded = NULL
) {
  shiny::moduleServer(id, function(input, output, session) {

    # ---- Resolve user-scoped paths ----
    # The caller (server.R) gates this module's construction on auth$logged,
    # so auth$user_dir / auth$email hold their finalized values by the time
    # we run. We snapshot them once at init — within a session they don't
    # change, and the SessionStore + downstream controllers expect plain
    # character paths, not reactives.
    #
    # Copilot chats and uploaded docs are personal, so they always live in
    # the per-user folder <pgx_dir>/<email> when an email is known —
    # independent of ENABLE_USERDIR. That flag only governs *dataset*
    # storage; when it is off it collapses auth$user_dir to the shared
    # PGX.DIR, which would otherwise scatter every user's chats/docs into one
    # place. Deriving the path from <pgx_dir>/<email> keeps Copilot pinned to
    # the user's own folder either way (when ENABLE_USERDIR is on this is the
    # same path auth$user_dir already points at). With no email (anonymous)
    # we fall back to auth$user_dir.
    user_dir <- shiny::isolate(auth$user_dir)
    email    <- shiny::isolate(auth$email)
    if (is.null(user_dir) || !nzchar(user_dir)) {
      stop("CopilotBoardServer: auth$user_dir is not set", call. = FALSE)
    }
    if (!is.null(email) && nzchar(email) &&
        !is.null(pgx_dir) && nzchar(pgx_dir)) {
      user_dir <- file.path(pgx_dir, email)
    }
    chat_dir <- file.path(user_dir, "chats")
    docs_dir <- file.path(user_dir, "docs_sources")
    dir.create(chat_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(docs_dir, recursive = TRUE, showWarnings = FALSE)

    # --- Board
    OmicsBoard("board", pgx, title="AI Copilot", infotext = NULL)
    
    # ---- Dataset context card output ----
    output$dataset_info <- shiny::renderUI({
      shiny::validate(shiny::need(
        !is.null(pgx) && !is.null(pgx$X),
        "No dataset loaded yet."
      ))
      ##pgx <- local_pgx()

      n_samples  <- if (!is.null(pgx$X)) ncol(pgx$X) else "?"
      n_features <- if (!is.null(pgx$X)) nrow(pgx$X) else "?"
      organism  <- if (!is.null(pgx$organism)) pgx$organism else "unknown"
      dataname  <- if (!is.null(pgx$name)) pgx$name else "unknown"
      description  <- if (!is.null(pgx$description)) pgx$description else "none"            

      contrasts <- "none"
      if (!is.null(pgx$contrasts)) {
        ct <- colnames(pgx$contrasts)
        if (length(ct) > 3L) {
          contrasts <- paste(c(ct[1:3], paste0("+ ", length(ct) - 3L, " more")), collapse = ", ")
        } else {
          contrasts <- paste(ct, collapse = ", ")
        }
      }

      shiny::tags$div(
        style = "font-size: 0.85em; line-height: 1.6em; padding: 1.3em 0.6em;",
        shiny::tags$div(shiny::strong("Dataset: "),  dataname),
        shiny::tags$div(shiny::strong("Description: "),  description),
        shiny::tags$div(shiny::strong("Organism: "),  organism),                
        shiny::tags$div(shiny::strong("Samples: "),   n_samples),
        shiny::tags$div(shiny::strong(tspan("Genes: ")), n_features),
        shiny::tags$div(shiny::strong("Contrasts: "), contrasts)
      )
    })

    # ---- Reactive state (board-owned) ----
    agent_rv                  <- shiny::reactiveVal(NULL)
    tier                      <- shiny::reactiveVal(tiers[[1]])
    style                     <- shiny::reactiveVal(COPILOT_STYLES[[1]])
    custom                    <- shiny::reactiveVal("")
    run_status                <- shiny::reactiveVal("idle")
    history_invalidation_tick <- shiny::reactiveVal(0L)
    restore_inflight          <- shiny::reactiveVal(NULL)
    pgx_loaded_event          <- shiny::reactiveVal(NULL)
    chat_event_rv             <- shiny::reactiveVal(NULL)

    # ---- Non-reactive state ----
    # Chats live in a per-user db so enterprise members sharing a chats/
    # folder don't collide. 
    chat_db_name <- .ai_tel_sessions_db_name(email)
    .detect_legacy_sessions(chat_dir, email)
    chat_store <- omicsagentovi::SessionStore(
      session_dir = chat_dir,
      db_name     = chat_db_name
    )

    # ---- Evidence module (constructed first so $append_artifact is available
    # when run/restore bindings factories run) ----
    evidence <- CopilotEvidenceServer("evidence", local_pgx = shiny::reactive(pgx))
    reports  <- CopilotReportsServer("reports", pgx = pgx)

    .set_agent_tools_enabled <- function(agent, enabled) {
      if (is.null(agent)) return(agent)
      tryCatch({
        agent@runtime$tools_enabled <- isTRUE(enabled)
        agent
      }, error = function(e) {
        log_info("copilot.tools.toggle_failed", msg = conditionMessage(e))
        agent
      })
    }

    # ---- Save controller ----
    save_ctrl <- copilot_save_controller(
      store                     = chat_store,
      agent                     = agent_rv,
      history_invalidation_tick = history_invalidation_tick,
      session                   = session,
      style                     = style,
      custom                    = custom
    )

    # ---- Telemetry: collect on save tick ----
    # The save controller bumps history_invalidation_tick AFTER session_save()
    # commits (post onFlushed), so by the time this fires the assistant turns
    # are on disk and a read-only scan sees them. The collect call is a tiny
    # SYNCHRONOUS main-thread op — deliberately NO future_promise / ExtendedTask:
    # the SQLite pointer + Shiny session closures would break under future
    # serialization (bad_weak_ptr), exactly as copilot_save_controller documents.
    # tryCatch keeps any telemetry failure from disrupting the chat UI.
    shiny::observeEvent(history_invalidation_tick(), {
      if (is.null(email) || length(email) != 1L || is.na(email) ||
          !nzchar(email)) {
        return(invisible(NULL))
      }
      tryCatch(
        ai_telemetry_collect_chat(session_dir = chat_dir, user_email = email),
        error = function(e) {
          log_info("copilot.telemetry_collect_failed", msg = conditionMessage(e))
        }
      )
    }, ignoreInit = TRUE)

    # ---- Restore controller ----
    restore_ctrl <- copilot_restore_controller(
      store            = chat_store,
      agent            = agent_rv,
      run_status       = run_status,
      restore_inflight = restore_inflight,
      bindings_factory = function() build_run_bindings(
                           session          = session,
                           evidence_api     = list(append_artifact = evidence$append_artifact),
                           docs_dir         = docs_dir,
                           data_dir         = pgx_dir,
                           pgx_loaded_event = pgx_loaded_event
                         ),
      local_pgx        = shiny::reactive(pgx),
      data_dir         = pgx_dir,
      evidence         = evidence,
      chat_event       = chat_event_rv,
      session          = session
    )

    # ---- Docs panel (constructed before run_ctrl so doc_context can be wired) ----
    docs <- CopilotDocsServer("docs", docs_dir = docs_dir)

    # ---- Chat module ----
    tier_choices_rx <- shiny::reactive({
      stats::setNames(tiers, vapply(tiers, copilot_tier_label, character(1)))
    })
    style_choices_rx <- shiny::reactive(COPILOT_STYLE_LABELS)
    chat_mod <- CopilotChatServer(
      "chat",
      on_user_message = function(text) {
        run_ctrl$dispatch(run_request_ask(text, show_user_msg = TRUE))
      },
      chat_event    = shiny::reactive(chat_event_rv()),
      run_status    = shiny::reactive(run_status()),
      on_abort      = function(reason) {
        run_ctrl$dispatch(run_request_abort(reason))
      },
      tier_choices  = tier_choices_rx,
      current_tier  = tier,
      style_choices = style_choices_rx,
      current_style = style,
      starters      = shiny::reactive(COPILOT_STARTERS)
    )

    # ---- Run controller ----
    run_ctrl <- copilot_run_controller(
      agent                = agent_rv,
      run_status           = run_status,
      tier                 = tier,
      save_ctrl            = save_ctrl,
      restore_ctrl         = restore_ctrl,
      evidence             = evidence,
      store                = chat_store,
      chat_event           = chat_event_rv,
      chat_on_tool_request = chat_mod$on_tool_request,
      pgx                  = pgx,
      pgx_dir              = pgx_dir,
      docs_dir             = docs_dir,
      pgx_loaded_event     = pgx_loaded_event,
      maxturns             = maxturns,
      session              = session,
      style                = style,
      custom               = custom,
      report_context       = reports,
      doc_context          = docs,
      tools_enabled        = reports$tools_enabled
    )

    # ---- Panel modules (datasets / history) ----
    # docs is constructed above so doc_context can be passed into run_ctrl.
    datasets <- CopilotDatasetsServer("datasets", pgx_dir = pgx_dir)
    history  <- CopilotHistoryServer(
      "history",
      session_dir               = chat_dir,
      db_name                   = chat_db_name,
      history_invalidation_tick = shiny::reactive(history_invalidation_tick())
    )

    shiny::observeEvent(pgx$name, {
      reports$refresh()
    }, ignoreNULL = FALSE)

    shiny::observeEvent(pgx$ai, {
      reports$refresh()
    }, ignoreNULL = FALSE)

    shiny::observeEvent(reports$tools_enabled(), {
      current <- shiny::isolate(agent_rv())
      if (!is.null(current)) {
        agent_rv(.set_agent_tools_enabled(current, reports$tools_enabled()))
      }
    }, ignoreInit = FALSE)

    # ---- Run dispatch observers (new chat / tier change) ----
    shiny::observeEvent(input$new_chat, {
      run_ctrl$dispatch(run_request_new_chat())
    }, ignoreInit = TRUE)
    # Tier changes come from the chat module's popover radioButtons (living in
    # the "chat" namespace) via the returned tier_clicked reactiveVal — the old
    # observeEvent(input$tier) was reading the board namespace and never fired.
    shiny::observeEvent(chat_mod$tier_clicked(), {
      new_tier <- chat_mod$tier_clicked()
      shiny::req(new_tier)
      if (identical(new_tier, shiny::isolate(tier()))) return()
      run_ctrl$dispatch(run_request_tier(new_tier))
    }, ignoreInit = TRUE)

    # Style + custom-text observers: rebuild the system prompt from fragments
    # and live-swap it onto the active ellmer Chat (via $set_system_prompt)
    # so the change takes effect on the next user message — no full reset
    # required (chat history is preserved). The S7 `system_prompt` slot is
    # updated in lockstep so the next session_save() persists the new prompt.
    .apply_style_to_live_agent <- function(agent, new_style, new_custom) {
      if (is.null(agent)) return(agent)
      new_prompt <- tryCatch(
        omicsagentovi::ovi_build_system_prompt(style = new_style, custom = new_custom),
        error = function(e) {
          log_info("copilot.style.build_failed", msg = conditionMessage(e))
          NULL
        }
      )
      if (is.null(new_prompt)) return(agent)
      tryCatch(agent@chat$set_system_prompt(new_prompt),
               error = function(e) {
                 log_info("copilot.style.live_swap_failed", msg = conditionMessage(e))
               })
      agent <- tryCatch(S7::set_props(agent, system_prompt = new_prompt),
                        error = function(e) agent)
      tryCatch(omicsagentovi::session_mark_changed(agent, source = "style"),
               error = function(e) agent)
    }

    shiny::observeEvent(chat_mod$style_clicked(), {
      new_style <- chat_mod$style_clicked()
      shiny::req(new_style)
      if (identical(new_style, shiny::isolate(style()))) return()
      style(new_style)
      current <- shiny::isolate(agent_rv())
      if (!is.null(current)) {
        agent_rv(.apply_style_to_live_agent(
          current, new_style, shiny::isolate(custom())
        ))
      }
    }, ignoreInit = TRUE)

    shiny::observeEvent(chat_mod$custom_text(), {
      new_custom <- chat_mod$custom_text() %||% ""
      if (identical(new_custom, shiny::isolate(custom()))) return()
      custom(new_custom)
      current <- shiny::isolate(agent_rv())
      if (!is.null(current)) {
        agent_rv(.apply_style_to_live_agent(
          current, shiny::isolate(style()), new_custom
        ))
      }
    }, ignoreInit = TRUE, ignoreNULL = FALSE)

    # ---- Dataset-change observers (route through run_ctrl$apply_dataset) ----
    # Source 1: pgx_loaded_event trampoline (manage_pgx tool via notification_sink)
    shiny::observeEvent(pgx_loaded_event(), {
      shiny::req(!is.null(pgx_loaded_event()))
      p <- pgx_loaded_event()
      # Use the bridge's resolved absolute path (omicsagentovi >= 0.4.1
      # adds dataset_path to the notification_sink payload). p$name_arg
      # is the raw LLM input — could be a basename — and must not be
      # stored as a locator path or restore will fail to re-load.
      # apply_dataset() funnels pgx_val through copilot_normalize_pgx()
      # internally — no need to snapshot here.
      run_ctrl$apply_dataset(
        pgx_val  = p$pgx,
        name     = p$dataset_name,
        path     = p$dataset_path,
        data_dir = p$data_dir
      )
      pgx_loaded_event(NULL)
    }, ignoreNULL = TRUE)

    # Source 2: global PGX from other boards. apply_dataset() funnels
    # pgx_val through copilot_normalize_pgx(), which handles the
    # reactiveValues -> classed list snapshot. We pass `pgx` directly
    # so there's only one snapshot site.
    shiny::observeEvent(pgx$name, {
      shiny::req(pgx$name, !is.null(pgx$X))
      if (!is.null(restore_inflight())) return()
      run_ctrl$apply_dataset(
        pgx_val  = pgx,
        name     = as.character(shiny::isolate(pgx$name)[[1]]),
        path     = NULL,
        data_dir = pgx_dir
      )
    }, ignoreNULL = TRUE)

    # ---- History panel: restore trigger ----
    shiny::observeEvent(history$on_restore(), {
      sid <- history$on_restore()
      shiny::req(sid)
      # Seed the style/custom reactives from the saved session metadata so
      # the radios reflect the restored session (and the next Agent build —
      # e.g. tier change after restore — uses the restored style).
      meta <- tryCatch(
        omicsagentovi::ovi_session_meta(session_id = sid, session_dir = chat_dir,
                                        db_name = chat_db_name),
        error = function(e) NULL
      )
      if (!is.null(meta)) {
        saved_style <- meta$style[[1L]]
        if (!is.na(saved_style) && nzchar(saved_style) &&
            saved_style %in% COPILOT_STYLES) {
          style(saved_style)
        }
        saved_custom <- meta$custom[[1L]]
        custom(if (is.na(saved_custom)) "" else saved_custom)
      }
      # Restoring a saved session starts a fresh chat from the user's POV —
      # re-arm both context tickboxes so previously-staged blocks can be
      # sent again under the restored agent.
      if (is.function(reports$reset_consumed)) reports$reset_consumed()
      if (is.function(docs$reset_consumed)) docs$reset_consumed()
      restore_ctrl$start(sid)
      history$on_restore(NULL)   # consume edge
    }, ignoreNULL = TRUE)

    # ---- History panel: delete trigger ----
    shiny::observeEvent(history$on_delete(), {
      sid <- history$on_delete()
      shiny::req(sid)
      tryCatch(
        omicsagentovi::ovi_session_delete(session_id = sid, session_dir = chat_dir,
                                          db_name = chat_db_name),
        error = function(e) {
          shiny::showNotification(
            paste("Delete failed:", conditionMessage(e)),
            type = "error", session = session
          )
        }
      )
      history_invalidation_tick(shiny::isolate(history_invalidation_tick()) + 1L)
      history$on_delete(NULL)
    }, ignoreNULL = TRUE)

    # ---- Greeting on first flush ----
    session$onFlushed(function() {
      chat_event_rv(list(
        type = "post", role = "assistant",
        text = copilot_msg("greeting")
      ))
    }, once = TRUE)

    invisible(NULL)
  })
}
