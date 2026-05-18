# copilot_server.R — Copilot Board Orchestrator
#
# Lifecycle owner for the Agent reactive. Wires module servers and lifecycle
# controllers. Does not implement run/save/restore logic — those live in their
# respective controllers.
#
# Phase 4 scope:
#   - Raw Agent reactive (agent_rv)
#   - Save + restore controllers wired (Phases 2 + 3)
#   - Chat module + run controller live (Phase 4)
#   - apply_dataset now lives in the run controller (moved from orchestrator)

# ---- Tier IDs (from contract_audit.md §1 / copilot-models.R:12-25) ----
.COPILOT_TIER_IDS <- c("copilot-default", "copilot-fast", "copilot-deep")

#' Copilot Board Orchestrator
#'
#' Wires module servers and controllers. Holds the agent reactive.
#' Does not contain run/save/restore logic — those live in their controllers.
#'
#' @param id         Module namespace id.
#' @param pgx        Global PGX reactiveValues shared across all boards.
#' @param pgx_dir    Path to directory containing .pgx dataset files.
#'                   Passed to RunBindings$data_dir and CopilotDatasetsServer.
#' @param chat_dir   Path for persisting chat sessions (SessionStore root).
#' @param docs_dir   Path for uploaded documents. Passed to RunBindings$docs_dir.
#' @param maxturns   Maximum user turns per session. Passed to run controller.
#' @param tiers      Character vector of tier identifiers (first = default).
#' @param is_data_loaded Optional reactiveVal; incremented after dataset load.
#'
#' @return NULL (side-effects only; called for Shiny registration).
#' @export
CopilotBoardServer <- function(
  id,
  pgx            = NULL,
  pgx_dir        = NULL,
  chat_dir,
  docs_dir,
  maxturns       = Inf,
  tiers          = "copilot-default",
  is_data_loaded = NULL
) {
  shiny::moduleServer(id, function(input, output, session) {

    # ---- Reactive state (board-owned) ----
    agent_rv                  <- shiny::reactiveVal(NULL)
    tier                      <- shiny::reactiveVal(tiers[[1]])
    run_status                <- shiny::reactiveVal("idle")
    history_invalidation_tick <- shiny::reactiveVal(0L)
    restore_inflight          <- shiny::reactiveVal(NULL)
    pgx_loaded_event          <- shiny::reactiveVal(NULL)
    chat_event_rv             <- shiny::reactiveVal(NULL)

    # ---- Non-reactive state ----
    chat_store <- omicsagentovi::SessionStore(session_dir = chat_dir)

    # ---- Phase 5: evidence module ----
    # Instantiated BEFORE run/restore controllers so $append_artifact is available
    # when bindings factories are called.
    evidence <- CopilotEvidenceServer("evidence", local_pgx = shiny::reactive(pgx))

    # ---- Phase 2: save controller ----
    save_ctrl <- copilot_save_controller(
      store                     = chat_store,
      agent                     = agent_rv,
      history_invalidation_tick = history_invalidation_tick,
      session                   = session
    )

    # ---- Phase 3: restore controller ----
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
    # TODO(phase 6): wire history$on_restore() -> restore_ctrl$start(sid)

    # ---- Phase 4: chat module ----
    tier_choices_rx <- shiny::reactive({
      # TODO(phase 6): use copilot_tier_label() for display labels.
      stats::setNames(tiers, tiers)
    })
    chat <- CopilotChatServer(
      "chat",
      on_user_message = function(text) {
        run_ctrl$dispatch(run_request_ask(text, show_user_msg = TRUE))
      },
      chat_event   = shiny::reactive(chat_event_rv()),
      run_status   = shiny::reactive(run_status()),
      tier_choices = tier_choices_rx
    )

    # ---- Phase 4: run controller ----
    run_ctrl <- copilot_run_controller(
      agent                = agent_rv,
      run_status           = run_status,
      tier                 = tier,
      save_ctrl            = save_ctrl,
      restore_ctrl         = restore_ctrl,
      evidence             = evidence,
      store                = chat_store,
      chat_event           = chat_event_rv,
      chat_on_tool_request = chat$on_tool_request,
      pgx                  = pgx,
      pgx_dir              = pgx_dir,
      docs_dir             = docs_dir,
      pgx_loaded_event     = pgx_loaded_event,
      maxturns             = maxturns,
      session              = session
    )

    # ---- Run dispatch observers (stop / new chat / tier change) ----
    shiny::observeEvent(input$stop_btn, {
      run_ctrl$dispatch(run_request_abort())
    })
    shiny::observeEvent(input$new_chat, {
      run_ctrl$dispatch(run_request_new_chat())
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$tier, {
      shiny::req(input$tier)
      if (identical(input$tier, shiny::isolate(tier()))) return()
      run_ctrl$dispatch(run_request_tier(input$tier))
    }, ignoreInit = TRUE)

    # ---- Dataset-change observers (route through run_ctrl$apply_dataset) ----
    # Source 1: pgx_loaded_event trampoline (manage_pgx tool via notification_sink)
    shiny::observeEvent(pgx_loaded_event(), {
      shiny::req(!is.null(pgx_loaded_event()))
      p <- pgx_loaded_event()
      run_ctrl$apply_dataset(
        pgx_val  = p$pgx,
        name     = p$dataset_name,
        path     = p$name_arg,
        data_dir = p$data_dir
      )
      pgx_loaded_event(NULL)
    }, ignoreNULL = TRUE)

    # Source 2: global PGX from other boards
    # Materialise the reactiveValues into a plain list with class "pgx" before
    # handing it to the agent — tools run outside the reactive domain and the
    # package's downstream code (.ovi_pgx_name + omicspgxmcp plot tools) expects
    # a plain list, not a reactiveValues handle.
    shiny::observeEvent(pgx$name, {
      shiny::req(pgx$name, !is.null(pgx$X))
      if (!is.null(restore_inflight())) return()
      snapshot <- shiny::reactiveValuesToList(pgx)
      class(snapshot) <- unique(c("pgx", class(snapshot)))
      run_ctrl$apply_dataset(
        pgx_val  = snapshot,
        name     = as.character(pgx$name[[1]]),
        path     = NULL,
        data_dir = pgx_dir
      )
    }, ignoreNULL = TRUE)

    # TODO(phase 6): Source 3 — datasets panel selection.

    # ---- Greeting on first flush ----
    session$onFlushed(function() {
      # TODO(phase 6): copilot_msg("greeting")
      chat_event_rv(list(
        type = "post", role = "assistant",
        text = "Hi — load a dataset and ask me anything about your experiment."
      ))
    }, once = TRUE)

    # ---- Phase 6: datasets/history/docs modules ----
    # datasets <- CopilotDatasetsServer("datasets", ...)
    # history  <- CopilotHistoryServer("history",  ...)
    # docs     <- CopilotDocsServer("docs",         ...)

    invisible(NULL)
  })
}
