## test-copilot-restore-controller.R
##
## Tests for copilot_restore_controller() and its internal helpers.
## Phase 4 refactor: shinychat I/O replaced by chat_event reactiveVal writes.
##
## Uses local_mocked_bindings() to stub omicsagentovi package calls.

# ---- Load omicsagentovi from source if the installed version is too old ----
.ovi_source_dir <- Sys.getenv("OVI_SOURCE_DIR", "/Users/santiago/projects/BigOmics/omicsagentovi")
if (packageVersion("omicsagentovi") < "0.4.0" && dir.exists(.ovi_source_dir)) {
  pkgload::load_all(.ovi_source_dir, quiet = TRUE)
}

if (packageVersion("omicsagentovi") < "0.4.0") {
  skip("omicsagentovi >= 0.4.0 required (agent_refactor branch). Set OVI_SOURCE_DIR.")
}

# ---- Source dependencies ----
# Historical name was `board.copilot`; the live module is `app_copilot`.
# Prefer the current name and fall back to the legacy one so test fixtures
# from older branches keep working.
.board_dir <- if (dir.exists("components/app_copilot/R")) {
  "components/app_copilot/R"
} else if (dir.exists("../../components/app_copilot/R")) {
  "../../components/app_copilot/R"
} else if (dir.exists("components/board.copilot/R")) {
  "components/board.copilot/R"
} else {
  "../../components/board.copilot/R"
}

source(file.path(.board_dir, "copilot_options.R"),  local = TRUE)
source(file.path(.board_dir, "copilot_messages.R"), local = TRUE)
source(file.path(.board_dir, "copilot_logger.R"),   local = TRUE)
source(file.path(.board_dir, "copilot_bindings.R"),           local = TRUE)
# Centralised pgx normaliser; the restore controller funnels every
# leak-prone pgx source through it before handing off to agent_set_pgx.
source(file.path(.board_dir, "copilot_pgx_normalize.R"),       local = TRUE)
# Context-block staging providers (current_dataset, future AI report, etc.)
# are invoked by .complete_restore after agent reconstruction.
source(file.path(.board_dir, "copilot_context_blocks.R"),      local = TRUE)
source(file.path(.board_dir, "copilot_restore_controller.R"), local = TRUE)

library(omicsagentovi)

# ---- Shared helpers ----

.make_stub_chat <- function() {
  list(
    get_turns  = function() list(),
    set_turns  = function(x) invisible(NULL)
  )
}

make_stub_agent <- function() {
  omicsagentovi::Agent(tools = list(), chat = .make_stub_chat())
}

make_fake_store <- function(session_dir = tempdir()) {
  omicsagentovi::SessionStore(session_dir = session_dir)
}

make_transcript_record <- function(role = "user", content = "Hello", visible = TRUE) {
  omicsagentovi::TranscriptRecord(
    idx          = 1L,
    role         = role,
    content_text = content,
    visibility   = if (visible) "visible" else "hidden"
  )
}

# Event-capturing chat_event writer for testing.
# Returns the underlying reactiveVal and an accessor for the events pushed.
make_chat_event_capture <- function() {
  events <- list()
  rv <- shiny::reactiveVal(NULL)
  writer <- function(ev) {
    if (missing(ev)) return(rv())
    events[[length(events) + 1L]] <<- ev
    rv(ev)
  }
  class(writer) <- c("reactiveVal", class(writer))
  list(rv = rv, writer = writer, events = function() events)
}

# ===========================================================================
# Unit: .replay_transcript — pushes one replay event with filtered records
# ===========================================================================

test_that(".replay_transcript with 0 records still emits one replay event with clear_first=TRUE", {
  # Restore-controller invariant: every call to .replay_transcript MUST
  # emit exactly one replay event, even with zero records, because the
  # event is what clears the "Restoring…" placeholder from chat. A
  # separate `clear` event would coalesce with this one at the
  # reactiveVal layer and be lost — see the "append on restore" bug.
  agent <- make_stub_agent()
  cap <- make_chat_event_capture()

  local_mocked_bindings(
    session_transcript = function(agent_session, view = "user") list(),
    .package = "omicsagentovi"
  )

  n <- .replay_transcript(agent, cap$writer)
  expect_equal(n, 0L)
  expect_length(cap$events(), 1L)
  ev <- cap$events()[[1]]
  expect_equal(ev$type, "replay")
  expect_true(isTRUE(ev$clear_first))
  expect_length(ev$records, 0L)
})

test_that(".replay_transcript with N visible records pushes a single replay event with those records", {
  agent <- make_stub_agent()
  cap <- make_chat_event_capture()

  records <- list(
    make_transcript_record("user",      "First user message",     TRUE),
    make_transcript_record("assistant", "First assistant reply",  TRUE),
    make_transcript_record("user",      "Second user message",    TRUE)
  )

  local_mocked_bindings(
    session_transcript = function(agent_session, view = "user") records,
    .package = "omicsagentovi"
  )

  n <- .replay_transcript(agent, cap$writer)

  expect_equal(n, 3L)
  expect_length(cap$events(), 1L)

  ev <- cap$events()[[1]]
  expect_equal(ev$type, "replay")
  expect_true(isTRUE(ev$clear_first))
  expect_length(ev$records, 3L)
  expect_equal(vapply(ev$records, function(r) r@role, character(1)),
               c("user", "assistant", "user"))
})

test_that(".replay_transcript drops records whose content_text is empty", {
  agent <- make_stub_agent()
  cap <- make_chat_event_capture()

  records <- list(
    make_transcript_record("user",      "Non-empty message",  TRUE),
    make_transcript_record("assistant", "",                   TRUE),  # dropped
    make_transcript_record("user",      "Another message",    TRUE)
  )

  local_mocked_bindings(
    session_transcript = function(agent_session, view = "user") records,
    .package = "omicsagentovi"
  )

  n <- .replay_transcript(agent, cap$writer)

  expect_equal(n, 2L)
  ev <- cap$events()[[1]]
  expect_length(ev$records, 2L)
})

# ===========================================================================
# Failure handler — chat_event receives a `post` event with the failure text
# ===========================================================================

test_that(".handle_restore_failure on model-fatal error: restore_inflight cleared and status failed", {
  store            <- make_fake_store()
  agent_rv         <- shiny::reactiveVal(make_stub_agent())
  run_status_rv    <- shiny::reactiveVal("idle")
  restore_inflight <- shiny::reactiveVal(NULL)
  chat_event_rv    <- shiny::reactiveVal(NULL)

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_restore_controller(
        store            = store,
        agent            = agent_rv,
        run_status       = run_status_rv,
        restore_inflight = restore_inflight,
        bindings_factory = function() {
          stop("ovi_restore: cannot restore session 'x' — model 'gpt-99' is not resolvable")
        },
        local_pgx        = shiny::reactive(NULL),
        data_dir         = tempdir(),
        evidence         = NULL,
        chat_event       = chat_event_rv,
        session          = session
      )
    },
    expr = {
      ctrl$start("test-session-id")
      session$flushReact()

      expect_null(shiny::isolate(restore_inflight()))
      expect_equal(shiny::isolate(ctrl$status()), "failed")

      ev <- shiny::isolate(chat_event_rv())
      expect_equal(ev$type, "post")
      expect_match(ev$text, "Restore failed", fixed = TRUE)
    }
  )
})

test_that(".handle_restore_failure on no-saved-session error: chat_event posts failure hint", {
  store            <- make_fake_store()
  agent_rv         <- shiny::reactiveVal(make_stub_agent())
  run_status_rv    <- shiny::reactiveVal("idle")
  restore_inflight <- shiny::reactiveVal(NULL)
  chat_event_rv    <- shiny::reactiveVal(NULL)

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_restore_controller(
        store            = store,
        agent            = agent_rv,
        run_status       = run_status_rv,
        restore_inflight = restore_inflight,
        bindings_factory = function() {
          stop("ovi_restore: no saved session with id 'missing-id'.")
        },
        local_pgx        = shiny::reactive(NULL),
        data_dir         = tempdir(),
        evidence         = NULL,
        chat_event       = chat_event_rv,
        session          = session
      )
    },
    expr = {
      ctrl$start("missing-id")
      session$flushReact()

      expect_null(shiny::isolate(restore_inflight()))
      expect_equal(shiny::isolate(ctrl$status()), "failed")

      ev <- shiny::isolate(chat_event_rv())
      expect_equal(ev$type, "post")
      expect_match(ev$text, "Restore failed", fixed = TRUE)
    }
  )
})

# ===========================================================================
# start() guards
# ===========================================================================

test_that("start() while run_status is 'streaming' does NOT invoke task", {
  store            <- make_fake_store()
  agent_rv         <- shiny::reactiveVal(make_stub_agent())
  run_status_rv    <- shiny::reactiveVal("streaming")
  restore_inflight <- shiny::reactiveVal(NULL)
  chat_event_rv    <- shiny::reactiveVal(NULL)

  ovi_restore_called <- FALSE
  local_mocked_bindings(
    ovi_restore = function(...) {
      ovi_restore_called <<- TRUE
      stop("should not be called")
    },
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_restore_controller(
        store            = store,
        agent            = agent_rv,
        run_status       = run_status_rv,
        restore_inflight = restore_inflight,
        bindings_factory = function() build_run_bindings(
                             session          = session,
                             evidence_api     = NULL,
                             pgx_loaded_event = NULL
                           ),
        local_pgx        = shiny::reactive(NULL),
        data_dir         = tempdir(),
        evidence         = NULL,
        chat_event       = chat_event_rv,
        session          = session
      )
    },
    expr = {
      ctrl$start("some-session-id")
      session$flushReact()

      expect_null(shiny::isolate(restore_inflight()))
      expect_false(ovi_restore_called)
      expect_equal(shiny::isolate(ctrl$status()), "idle")
    }
  )
})

test_that("start() while another restore is in flight returns silently", {
  store            <- make_fake_store()
  agent_rv         <- shiny::reactiveVal(make_stub_agent())
  run_status_rv    <- shiny::reactiveVal("idle")
  restore_inflight <- shiny::reactiveVal("existing-session-id")
  chat_event_rv    <- shiny::reactiveVal(NULL)

  ovi_restore_called <- FALSE
  local_mocked_bindings(
    ovi_restore = function(...) {
      ovi_restore_called <<- TRUE
      stop("should not be called")
    },
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_restore_controller(
        store            = store,
        agent            = agent_rv,
        run_status       = run_status_rv,
        restore_inflight = restore_inflight,
        bindings_factory = function() build_run_bindings(
                             session          = session,
                             evidence_api     = NULL,
                             pgx_loaded_event = NULL
                           ),
        local_pgx        = shiny::reactive(NULL),
        data_dir         = tempdir(),
        evidence         = NULL,
        chat_event       = chat_event_rv,
        session          = session
      )
    },
    expr = {
      ctrl$start("new-session-id")
      session$flushReact()

      expect_equal(shiny::isolate(restore_inflight()), "existing-session-id")
      expect_false(ovi_restore_called)
    }
  )
})

test_that("start() happy path: chat clear + post events pushed, evidence cleared, task invoked", {
  store            <- make_fake_store()
  agent_rv         <- shiny::reactiveVal(make_stub_agent())
  run_status_rv    <- shiny::reactiveVal("idle")
  restore_inflight <- shiny::reactiveVal(NULL)
  chat_event_rv    <- shiny::reactiveVal(NULL)

  events <- list()
  push_writer <- function(ev) {
    events[[length(events) + 1L]] <<- ev
    chat_event_rv(ev)
  }

  evidence_cleared   <- FALSE
  ovi_restore_called <- FALSE
  fake_evidence <- list(clear = function() { evidence_cleared <<- TRUE })

  local_mocked_bindings(
    ovi_restore = function(session_id, session_dir, bindings, restore_pgx) {
      ovi_restore_called <<- TRUE
      make_stub_agent()
    },
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_restore_controller(
        store            = store,
        agent            = agent_rv,
        run_status       = run_status_rv,
        restore_inflight = restore_inflight,
        bindings_factory = function() build_run_bindings(
                             session          = session,
                             evidence_api     = NULL,
                             pgx_loaded_event = NULL
                           ),
        local_pgx        = shiny::reactive(NULL),
        data_dir         = tempdir(),
        evidence         = fake_evidence,
        chat_event       = push_writer,
        session          = session
      )
    },
    expr = {
      ctrl$start("test-session-id")
      session$flushReact()

      expect_equal(shiny::isolate(restore_inflight()), "test-session-id")
      expect_true(evidence_cleared)

      # restore_controller emits a single atomic `reset` event (clear +
      # placeholder post in one observer firing). Two back-to-back writes
      # to chat_event would coalesce at the reactiveVal layer and the
      # `clear` would be lost — that was the "append on restore" bug.
      expect_true(length(events) >= 1L)
      expect_equal(events[[1]]$type, "reset")
      expect_true(nzchar(events[[1]]$text %||% ""))

      expect_equal(shiny::isolate(ctrl$status()), "restoring")
    }
  )
})

# ===========================================================================
# Result observer paths (pre-resolution invariants only — ExtendedTask
# does not run real futures in testServer)
# ===========================================================================

test_that("after task starts with NULL local_pgx: inflight set, status restoring", {
  store            <- make_fake_store()
  agent_rv         <- shiny::reactiveVal(make_stub_agent())
  run_status_rv    <- shiny::reactiveVal("idle")
  restore_inflight <- shiny::reactiveVal(NULL)
  chat_event_rv    <- shiny::reactiveVal(NULL)

  agent_set_pgx_called <- FALSE

  local_mocked_bindings(
    ovi_restore = function(session_id, session_dir, bindings, restore_pgx) make_stub_agent(),
    agent_set_pgx = function(agent, pgx, ...) {
      agent_set_pgx_called <<- TRUE
      agent
    },
    session_transcript = function(agent_session, view = "user") list(),
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_restore_controller(
        store            = store,
        agent            = agent_rv,
        run_status       = run_status_rv,
        restore_inflight = restore_inflight,
        bindings_factory = function() build_run_bindings(
                             session          = session,
                             evidence_api     = NULL,
                             pgx_loaded_event = NULL
                           ),
        local_pgx        = shiny::reactive(NULL),
        data_dir         = tempdir(),
        evidence         = NULL,
        chat_event       = chat_event_rv,
        session          = session
      )
    },
    expr = {
      ctrl$start("test-session-id")
      session$flushReact()
      expect_equal(shiny::isolate(restore_inflight()), "test-session-id")
      expect_false(agent_set_pgx_called)
    }
  )
})

test_that("non-NULL local_pgx triggers the sync fast path: restore completes inside session$onFlushed", {
  # Sync fast path contract: when local_pgx is non-NULL, ovi_restore runs
  # inside session$onFlushed and .complete_restore fires immediately after.
  # By the time flushReact() returns, the controller has already cleared
  # restore_inflight and dropped status back to "idle". (The pre-fix
  # assertion checked the "restoring" mid-state — that only holds in the
  # async path; for the sync path the right invariant is post-completion.)
  store            <- make_fake_store()
  restored_agent   <- make_stub_agent()
  agent_rv         <- shiny::reactiveVal(make_stub_agent())
  run_status_rv    <- shiny::reactiveVal("idle")
  restore_inflight <- shiny::reactiveVal(NULL)
  chat_event_rv    <- shiny::reactiveVal(NULL)
  # copilot_normalize_pgx now requires both $name and $X to consider a
  # value a usable PGX — incomplete shapes get treated as "no dataset
  # yet" and return NULL. Give the fake both fields so agent_set_pgx
  # actually fires on the sync path.
  fake_pgx <- list(name = "test-dataset", X = matrix(1, 2, 2))
  agent_set_pgx_called <- FALSE

  local_mocked_bindings(
    ovi_restore = function(session_id, session_dir, bindings, restore_pgx) restored_agent,
    agent_set_pgx = function(agent, pgx, ...) {
      agent_set_pgx_called <<- TRUE
      agent
    },
    session_transcript = function(agent_session, view = "user") list(),
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_restore_controller(
        store            = store,
        agent            = agent_rv,
        run_status       = run_status_rv,
        restore_inflight = restore_inflight,
        bindings_factory = function() build_run_bindings(
                             session          = session,
                             evidence_api     = NULL,
                             pgx_loaded_event = NULL
                           ),
        local_pgx        = shiny::reactive(fake_pgx),
        data_dir         = tempdir(),
        evidence         = NULL,
        chat_event       = chat_event_rv,
        session          = session
      )
    },
    expr = {
      ctrl$start("test-session-id")
      session$flushReact()
      # Sync path complete: inflight cleared, status idle, pgx injected,
      # restored agent committed.
      expect_null(shiny::isolate(restore_inflight()))
      expect_equal(shiny::isolate(ctrl$status()), "idle")
      expect_true(agent_set_pgx_called)
      expect_identical(shiny::isolate(agent_rv()), restored_agent)
    }
  )
})

# ===========================================================================
# Failure invariants
# ===========================================================================

test_that("failure during start() leaves agent() unchanged (previous_agent preserved)", {
  store            <- make_fake_store()
  original_agent   <- make_stub_agent()
  agent_rv         <- shiny::reactiveVal(original_agent)
  run_status_rv    <- shiny::reactiveVal("idle")
  restore_inflight <- shiny::reactiveVal(NULL)
  chat_event_rv    <- shiny::reactiveVal(NULL)

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_restore_controller(
        store            = store,
        agent            = agent_rv,
        run_status       = run_status_rv,
        restore_inflight = restore_inflight,
        bindings_factory = function() {
          stop("ovi_restore: cannot restore session 'x' — model 'old-model' is not resolvable")
        },
        local_pgx        = shiny::reactive(NULL),
        data_dir         = tempdir(),
        evidence         = NULL,
        chat_event       = chat_event_rv,
        session          = session
      )
    },
    expr = {
      ctrl$start("fatal-session-id")
      session$flushReact()

      expect_null(shiny::isolate(restore_inflight()))
      expect_equal(shiny::isolate(ctrl$status()), "failed")
      expect_identical(shiny::isolate(agent_rv()), original_agent)
    }
  )
})

test_that("agent_set_pgx throws during injection: restore still completes, agent committed without pgx", {
  # Resilience contract: if agent_set_pgx errors (e.g. omicspgx accessor
  # blows up on a malformed dataset), .complete_restore must log and
  # proceed — the rest of the restore (context block staging, transcript
  # replay) must still land, and the controller must reach idle. We
  # supply a well-shaped fake_pgx so the normaliser lets it through; the
  # MOCKED agent_set_pgx is what throws.
  store            <- make_fake_store()
  restored_agent   <- make_stub_agent()
  agent_rv         <- shiny::reactiveVal(make_stub_agent())
  run_status_rv    <- shiny::reactiveVal("idle")
  restore_inflight <- shiny::reactiveVal(NULL)
  chat_event_rv    <- shiny::reactiveVal(NULL)
  fake_pgx <- list(name = "bad-dataset", X = matrix(1, 2, 2))
  set_pgx_attempted <- FALSE

  local_mocked_bindings(
    ovi_restore = function(session_id, session_dir, bindings, restore_pgx) restored_agent,
    agent_set_pgx = function(agent, pgx, ...) {
      set_pgx_attempted <<- TRUE
      stop("agent_set_pgx: simulated failure")
    },
    session_transcript = function(agent_session, view = "user") list(),
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_restore_controller(
        store            = store,
        agent            = agent_rv,
        run_status       = run_status_rv,
        restore_inflight = restore_inflight,
        bindings_factory = function() build_run_bindings(
                             session          = session,
                             evidence_api     = NULL,
                             pgx_loaded_event = NULL
                           ),
        local_pgx        = shiny::reactive(fake_pgx),
        data_dir         = tempdir(),
        evidence         = NULL,
        chat_event       = chat_event_rv,
        session          = session
      )
    },
    expr = {
      ctrl$start("test-session-id")
      session$flushReact()

      # Injection was attempted (the normaliser accepted the fake_pgx).
      expect_true(set_pgx_attempted)
      # And restoration still completed cleanly — failure of the pgx
      # inject must not strand inflight or block agent commit.
      expect_null(shiny::isolate(restore_inflight()))
      expect_equal(shiny::isolate(ctrl$status()), "idle")
      expect_identical(shiny::isolate(agent_rv()), restored_agent)
    }
  )
})

# ===========================================================================
# Public surface shape
# ===========================================================================

test_that("copilot_restore_controller returns start function and status reactive", {
  store            <- make_fake_store()
  agent_rv         <- shiny::reactiveVal(NULL)
  run_status_rv    <- shiny::reactiveVal("idle")
  restore_inflight <- shiny::reactiveVal(NULL)
  chat_event_rv    <- shiny::reactiveVal(NULL)

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_restore_controller(
        store            = store,
        agent            = agent_rv,
        run_status       = run_status_rv,
        restore_inflight = restore_inflight,
        bindings_factory = function() build_run_bindings(
                             session          = session,
                             evidence_api     = NULL,
                             pgx_loaded_event = NULL
                           ),
        local_pgx        = shiny::reactive(NULL),
        data_dir         = tempdir(),
        evidence         = NULL,
        chat_event       = chat_event_rv,
        session          = session
      )
    },
    expr = {
      expect_true(is.function(ctrl$start))
      expect_true(shiny::is.reactive(ctrl$status))
      expect_equal(shiny::isolate(ctrl$status()), "idle")
    }
  )
})
