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
#' @param chat_dir   Path for persisting chat sessions (SessionStore root).
#' @param docs_dir   Path for uploaded documents.
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
  chat_dir,
  docs_dir,
  maxturns       = Inf,
  tiers          = COPILOT_TIERS,
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

    # ---- Evidence module (constructed first so $append_artifact is available
    # when run/restore bindings factories run) ----
    evidence <- CopilotEvidenceServer("evidence", local_pgx = shiny::reactive(pgx))

    # ---- Save controller ----
    save_ctrl <- copilot_save_controller(
      store                     = chat_store,
      agent                     = agent_rv,
      history_invalidation_tick = history_invalidation_tick,
      session                   = session
    )

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

    # ---- Chat module ----
    tier_choices_rx <- shiny::reactive({
      stats::setNames(tiers, vapply(tiers, copilot_tier_label, character(1)))
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
      chat_on_tool_request = chat$on_tool_request,
      pgx                  = pgx,
      pgx_dir              = pgx_dir,
      docs_dir             = docs_dir,
      pgx_loaded_event     = pgx_loaded_event,
      maxturns             = maxturns,
      session              = session
    )

    # ---- Panel modules (datasets / history / docs) ----
    datasets <- CopilotDatasetsServer("datasets", pgx_dir = pgx_dir)
    history  <- CopilotHistoryServer(
      "history",
      session_dir               = chat_dir,
      history_invalidation_tick = shiny::reactive(history_invalidation_tick())
    )
    docs     <- CopilotDocsServer("docs", docs_dir = docs_dir)

    # ---- Run dispatch observers (new chat / tier change) ----
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

    # Source 2: global PGX from other boards. Materialise reactiveValues into
    # a plain list with class "pgx" — tools run outside the reactive domain
    # and downstream package code expects a plain list.
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

    # Source 3: human dataset picker (CopilotDatasetsServer selection)
    shiny::observeEvent(datasets$selected(), {
      path <- datasets$selected()
      shiny::req(path)
      if (!is.null(restore_inflight())) return()
      loaded <- tryCatch(
        {
          pgx_obj <- playbase::pgx.load(path)
          pgx_obj <- playbase::pgx.initialize(pgx_obj)
          pgx_obj
        },
        error = function(e) {
          shiny::showNotification(
            copilot_msg("switch_failed", msg = conditionMessage(e)),
            type = "error", session = session
          )
          NULL
        }
      )
      if (is.null(loaded)) return()
      # Push into global pgx so other boards see it.
      sync_rv_from_list(pgx, loaded)
      class(loaded) <- unique(c("pgx", class(loaded)))
      run_ctrl$apply_dataset(
        pgx_val  = loaded,
        name     = tools::file_path_sans_ext(basename(path)),
        path     = path,
        data_dir = pgx_dir
      )
      bigdash.showTabsGoToDataView(session)
      if (!is.null(is_data_loaded)) {
        is_data_loaded(shiny::isolate(is_data_loaded()) + 1L)
      }
    }, ignoreNULL = TRUE)

    # ---- History panel: restore trigger ----
    shiny::observeEvent(history$on_restore(), {
      sid <- history$on_restore()
      shiny::req(sid)
      restore_ctrl$start(sid)
      history$on_restore(NULL)   # consume edge
    }, ignoreNULL = TRUE)

    # ---- History panel: delete trigger ----
    shiny::observeEvent(history$on_delete(), {
      sid <- history$on_delete()
      shiny::req(sid)
      tryCatch(
        omicsagentovi::ovi_session_delete(session_id = sid, session_dir = chat_dir),
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
