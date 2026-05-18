# copilot_server.R — Copilot Board Orchestrator
#
# Lifecycle owner for the Agent reactive. Wires module servers and lifecycle
# controllers. Does not implement run/save/restore logic — those live in their
# respective controllers (landing in Phases 2–4).
#
# Phase 1 scope:
#   - Raw Agent reactive (agent_rv)
#   - pgx_loaded_event trampoline
#   - apply_dataset() single path for all three dataset-change sources
#   - Scaffolded stubs for Phases 2–6 (commented controllers + modules)

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
    agent_rv                  <- shiny::reactiveVal(NULL)  # raw omicsagentovi::Agent or NULL
    tier                      <- shiny::reactiveVal(tiers[[1]])
    run_status                <- shiny::reactiveVal("idle")  # "idle"|"streaming"|"settling"|"failed"
    history_invalidation_tick <- shiny::reactiveVal(0L)
    restore_inflight          <- shiny::reactiveVal(NULL)   # NULL or restoring session_id
    pgx_loaded_event          <- shiny::reactiveVal(NULL)   # trampoline for notification_sink

    # ---- Non-reactive state ----
    chat_store <- omicsagentovi::SessionStore(session_dir = chat_dir)

    # ---- apply_dataset() — single path for all three dataset-change sources ----
    # Sources:
    #   1. pgx_loaded_event trampoline (agent's manage_pgx tool via notification_sink)
    #   2. pgx$name reactive (global PGX from other boards)
    #   3. TODO(phase 6): datasets panel selection -> CopilotDatasetsServer
    #
    # Chat is preserved across switches (decision Q3 — overview.md §Decisions).
    apply_dataset <- function(payload) {
      if (is.null(agent_rv())) {
        # First dataset load: construct a fresh Agent.
        bindings <- build_run_bindings(
          session          = session,
          evidence_api     = NULL,  # TODO(phase 5): wire evidence$api here
          docs_dir         = docs_dir,
          data_dir         = pgx_dir,
          pgx_loaded_event = pgx_loaded_event
        )
        new_agent <- tryCatch(
          omicsagentovi::Agent(
            tier     = tier(),
            context  = omicsagentovi::RunContext(pgx = payload$pgx),
            bindings = bindings
          ),
          error = function(e) {
            shiny::showNotification(
              paste("Copilot: failed to create agent —", conditionMessage(e)),
              type = "error", session = session
            )
            NULL
          }
        )
        agent_rv(new_agent)
      } else {
        # Dataset switch: preserve Agent identity and chat history, update PGX.
        updated <- tryCatch(
          omicsagentovi::agent_set_pgx(
            agent_rv(),
            pgx          = payload$pgx,
            dataset_name = payload$dataset_name,
            dataset_path = payload$dataset_path,
            data_dir     = payload$data_dir
          ),
          error = function(e) {
            shiny::showNotification(
              paste("Copilot: failed to switch dataset —", conditionMessage(e)),
              type = "error", session = session
            )
            NULL
          }
        )
        if (!is.null(updated)) {
          agent_rv(updated)
        }
      }
    }

    # ---- Source 1: pgx_loaded_event trampoline ----
    # notification_sink writes payload here from the tool-execution thread.
    # We consume it on the Shiny thread and reset to NULL after.
    shiny::observeEvent(pgx_loaded_event(), {
      shiny::req(!is.null(pgx_loaded_event()))
      p <- pgx_loaded_event()
      apply_dataset(list(
        pgx          = p$pgx,
        dataset_name = p$dataset_name,
        dataset_path = p$name_arg,   # tool field name is name_arg; maps to dataset_path
        data_dir     = p$data_dir
      ))
      pgx_loaded_event(NULL)
    }, ignoreNULL = TRUE)

    # ---- Source 2: global PGX from other boards ----
    shiny::observeEvent(pgx$name, {
      shiny::req(pgx$name)
      if (!is.null(restore_inflight())) return()
      apply_dataset(list(
        pgx          = pgx,          # pass the reactiveValues handle; agent_set_pgx accepts it
        dataset_name = as.character(pgx$name[[1]]),
        dataset_path = NULL,
        data_dir     = pgx_dir
      ))
    }, ignoreNULL = TRUE)

    # TODO(phase 6): wire datasets$selected_dataset() -> apply_dataset()
    # When CopilotDatasetsServer is wired (Phase 6 rename):
    #   shiny::observeEvent(datasets$selected_dataset(), {
    #     shiny::req(!is.null(datasets$selected_dataset()))
    #     apply_dataset(list(
    #       pgx          = <loaded pgx from datasets server>,
    #       dataset_name = datasets$selected_dataset()$name,
    #       dataset_path = datasets$selected_dataset()$path,
    #       data_dir     = pgx_dir
    #     ))
    #   })

    # ---- Phase 2: save controller ----
    save_ctrl <- copilot_save_controller(
      store                     = chat_store,
      agent                     = agent_rv,
      history_invalidation_tick = history_invalidation_tick,
      session                   = session
    )
    # TODO(phase 4): replace with run_ctrl-emitted run_settled signal.
    # For now, save_ctrl$on_run_settled() is dormant — Phase 4 wires the trigger.

    # ---- Phase 3: restore controller ----
    # restore_ctrl <- copilot_restore_controller(
    #   bindings_factory = function() build_run_bindings(
    #                        session          = session,
    #                        evidence_api     = evidence$api,
    #                        docs_dir         = docs_dir,
    #                        data_dir         = pgx_dir,
    #                        pgx_loaded_event = pgx_loaded_event
    #                      ),
    #   store            = chat_store,
    #   agent_rv         = agent_rv,
    #   chat_api         = chat,
    #   evidence_api     = evidence,
    #   restore_inflight = restore_inflight,
    #   session          = session
    # )
    # shiny::observeEvent(history$on_restore(), {
    #   restore_ctrl$start(history$on_restore())
    # })

    # ---- Phase 4: run controller ----
    # run_ctrl <- copilot_run_controller(
    #   agent_rv         = agent_rv,
    #   store            = chat_store,
    #   save_ctrl        = save_ctrl,
    #   chat_api         = chat,
    #   evidence_api     = evidence,
    #   run_status       = run_status,
    #   session          = session
    # )
    # shiny::observeEvent(input$chat_user_input, {
    #   shiny::req(!identical(run_status(), "streaming"))
    #   run_ctrl$dispatch(list(kind = "ask", text = input$chat_user_input))
    # })
    # shiny::observeEvent(input$new_chat, {
    #   run_ctrl$dispatch(list(kind = "new_chat"))
    # }, ignoreInit = TRUE)
    # shiny::observeEvent(input$tier, {
    #   run_ctrl$dispatch(list(kind = "tier", new_tier = input$tier))
    # }, ignoreInit = TRUE)
    # shiny::observeEvent(input$stop_btn, {
    #   run_ctrl$dispatch(list(kind = "abort"))
    # })
    # shiny::observeEvent(input$ask_describe, {
    #   run_ctrl$dispatch(list(kind = "ask", text = .COPILOT_EXAMPLE_QUESTIONS$describe))
    # })
    # shiny::observeEvent(input$ask_findings, {
    #   run_ctrl$dispatch(list(kind = "ask", text = .COPILOT_EXAMPLE_QUESTIONS$findings))
    # })
    # shiny::observeEvent(input$ask_pathways, {
    #   run_ctrl$dispatch(list(kind = "ask", text = .COPILOT_EXAMPLE_QUESTIONS$pathways))
    # })
    # shiny::observeEvent(input$ask_biomarkers, {
    #   run_ctrl$dispatch(list(kind = "ask", text = .COPILOT_EXAMPLE_QUESTIONS$biomarkers))
    # })
    # shiny::observeEvent(input$ask_plot, {
    #   run_ctrl$dispatch(list(kind = "ask", text = .COPILOT_EXAMPLE_QUESTIONS$plot))
    # })

    # ---- Phase 4: chat module ----
    # chat <- CopilotChatServer("chat",
    #   agent_rv   = agent_rv,
    #   run_status = run_status,
    #   tier       = tier,
    #   tiers      = tiers
    # )

    # ---- Phase 5: evidence module ----
    # evidence <- CopilotEvidenceServer("evidence", pgx = pgx)
    # Re-wire build_run_bindings() after evidence is available:
    # each tier change or new_chat reconstructs bindings with evidence$api.

    # ---- Phase 6: datasets/history/docs modules ----
    # datasets <- CopilotDatasetsServer("datasets",
    #   pgx_dir        = pgx_dir,
    #   pgx            = pgx,
    #   is_data_loaded = is_data_loaded
    # )
    # history <- CopilotHistoryServer("history",
    #   chat_store       = chat_store,
    #   refresh_trigger  = history_invalidation_tick
    # )
    # docs <- CopilotDocsServer("docs", docs_dir = docs_dir)
    #
    # History delete observer:
    # shiny::observeEvent(history$on_delete(), {
    #   shiny::req(!is.null(history$on_delete()))
    #   omicsagentovi::ovi_session_delete(
    #     session_id   = history$on_delete(),
    #     session_dir  = chat_dir
    #   )
    #   history_invalidation_tick(history_invalidation_tick() + 1L)
    # })

    # ---- Tier selector initialisation ----
    session$onFlushed(function() {
      shiny::updateSelectInput(
        session = session,
        inputId = "tier",
        choices = setNames(tiers, tiers),  # TODO(phase 6): use copilot_tier_label()
        selected = tiers[[1]]
      )
    }, once = TRUE)

    invisible(NULL)
  })
}
