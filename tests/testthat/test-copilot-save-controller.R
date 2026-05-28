## test-copilot-save-controller.R
##
## Tests for copilot_save_controller() and its internal .prune_sessions helper.
##
## The controller now runs session_save() synchronously on the main thread,
## deferred via session$onFlushed(once = TRUE). shiny::testServer's
## session$flushReact() invokes onFlushed callbacks, so the save body fires
## inside the test harness without any future/ExtendedTask plumbing.
##
## Uses local_mocked_bindings() to stub omicsagentovi package calls.
## Does NOT require a live LLM or real SQLite store.

# ---- Load omicsagentovi from source if the installed version is too old ----
.ovi_source_dir <- Sys.getenv("OVI_SOURCE_DIR", "/Users/santiago/projects/BigOmics/omicsagentovi")
if (packageVersion("omicsagentovi") < "0.4.0" && dir.exists(.ovi_source_dir)) {
  pkgload::load_all(.ovi_source_dir, quiet = TRUE)
}

if (packageVersion("omicsagentovi") < "0.4.0") {
  skip("omicsagentovi >= 0.4.0 required (agent_refactor branch). Set OVI_SOURCE_DIR.")
}

# ---- Source the controller under test ----
.board_dir_shared <- if (dir.exists("components/board.copilot/R")) {
  "components/board.copilot/R"
} else {
  "../../components/board.copilot/R"
}
source(file.path(.board_dir_shared, "copilot_options.R"),  local = TRUE)
source(file.path(.board_dir_shared, "copilot_messages.R"), local = TRUE)
source(file.path(.board_dir_shared, "copilot_logger.R"),   local = TRUE)

.ctrl_path <- if (file.exists("components/board.copilot/R/copilot_save_controller.R")) {
  "components/board.copilot/R/copilot_save_controller.R"
} else {
  "../../components/board.copilot/R/copilot_save_controller.R"
}
source(.ctrl_path, local = TRUE)

library(omicsagentovi)

# ---- Helpers ----

.make_stub_chat <- function() {
  list(
    get_turns = function() list(),
    set_turns = function(x) invisible(NULL)
  )
}

make_stub_agent <- function() {
  omicsagentovi::Agent(tools = list(), chat = .make_stub_chat())
}

make_fake_store <- function(session_dir = tempdir()) {
  omicsagentovi::SessionStore(session_dir = session_dir)
}

# Build an Agent stand-in for session_save() to return — must carry a session
# with matching session_id so the result-write-back guard passes.
make_saved_agent_for <- function(original, bump = 1L) {
  saved_session <- S7::set_props(
    original@session,
    saved_generation = original@session@dirty_generation + bump
  )
  S7::set_props(original, session = saved_session)
}

# ===========================================================================
# Unit: .prune_sessions
# ===========================================================================

test_that(".prune_sessions with 5 sessions and max=3 deletes the 2 oldest", {
  store <- make_fake_store()

  fake_rows <- data.frame(
    session_id = c("s5", "s4", "s3", "s2", "s1"),  # newest -> oldest
    updated_at = as.numeric(c(5, 4, 3, 2, 1)),
    stringsAsFactors = FALSE
  )

  deleted <- character(0)

  local_mocked_bindings(
    ovi_sessions       = function(session_dir = NULL) fake_rows,
    ovi_session_delete = function(session_id, session_dir = NULL) {
      deleted <<- c(deleted, session_id)
      invisible(TRUE)
    },
    .package = "omicsagentovi"
  )

  .prune_sessions(store, max_history = 3L)
  expect_equal(sort(deleted), c("s1", "s2"))
})

test_that(".prune_sessions with sessions <= max_history does nothing", {
  store <- make_fake_store()

  fake_rows <- data.frame(
    session_id = c("s3", "s2", "s1"),
    updated_at = as.numeric(c(3, 2, 1)),
    stringsAsFactors = FALSE
  )

  deleted <- character(0)

  local_mocked_bindings(
    ovi_sessions       = function(session_dir = NULL) fake_rows,
    ovi_session_delete = function(session_id, session_dir = NULL) {
      deleted <<- c(deleted, session_id)
      invisible(TRUE)
    },
    .package = "omicsagentovi"
  )

  .prune_sessions(store, max_history = 3L)
  expect_equal(deleted, character(0))
})

test_that(".prune_sessions prune error is caught and does not propagate", {
  store <- make_fake_store()
  fake_rows <- data.frame(
    session_id = c("s2", "s1"),
    updated_at = as.numeric(c(2, 1)),
    stringsAsFactors = FALSE
  )
  local_mocked_bindings(
    ovi_sessions       = function(session_dir = NULL) fake_rows,
    ovi_session_delete = function(session_id, session_dir = NULL) stop("DB error"),
    .package = "omicsagentovi"
  )
  expect_no_error(.prune_sessions(store, max_history = 1L))
})

# ===========================================================================
# Controller construction
# ===========================================================================

test_that("copilot_save_controller returns on_run_settled function and status reactive", {
  store    <- make_fake_store()
  agent_rv <- shiny::reactiveVal(NULL)
  tick     <- shiny::reactiveVal(0L)

  shiny::testServer(
    function(input, output, session) {
      save_ctrl <<- copilot_save_controller(
        store                     = store,
        agent                     = agent_rv,
        history_invalidation_tick = tick,
        session                   = session
      )
    },
    expr = {
      expect_true(is.function(save_ctrl$on_run_settled))
      expect_true(shiny::is.reactive(save_ctrl$status))
      expect_equal(shiny::isolate(save_ctrl$status()), "idle")
    }
  )
})

# ===========================================================================
# on_run_settled: NULL agent -> no-op
# ===========================================================================

test_that("on_run_settled with NULL agent is a no-op", {
  store    <- make_fake_store()
  tick     <- shiny::reactiveVal(0L)
  agent_rv <- shiny::reactiveVal(NULL)

  shiny::testServer(
    function(input, output, session) {
      save_ctrl <<- copilot_save_controller(
        store                     = store,
        agent                     = agent_rv,
        history_invalidation_tick = tick,
        session                   = session
      )
    },
    expr = {
      expect_no_error(save_ctrl$on_run_settled())
      session$flushReact()
      expect_equal(shiny::isolate(save_ctrl$status()), "idle")
      expect_equal(shiny::isolate(tick()), 0L)
    }
  )
})

# ===========================================================================
# on_run_settled: clean agent -> no save
# ===========================================================================

test_that("on_run_settled with clean agent keeps status idle and tick unchanged", {
  store    <- make_fake_store()
  tick     <- shiny::reactiveVal(0L)
  agent    <- make_stub_agent()  # clean by construction
  agent_rv <- shiny::reactiveVal(agent)

  expect_false(omicsagentovi::session_is_dirty(agent@session))

  saved <- FALSE
  local_mocked_bindings(
    session_save = function(store, agent, ...) { saved <<- TRUE; agent },
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      save_ctrl <<- copilot_save_controller(
        store                     = store,
        agent                     = agent_rv,
        history_invalidation_tick = tick,
        session                   = session
      )
    },
    expr = {
      save_ctrl$on_run_settled()
      session$flushReact()
      expect_false(saved)
      expect_equal(shiny::isolate(save_ctrl$status()), "idle")
      expect_equal(shiny::isolate(tick()), 0L)
    }
  )
})

# ===========================================================================
# on_run_settled: dirty agent -> save runs on flush, status returns to idle,
# tick bumped, agent_rv updated with saved agent
# ===========================================================================

test_that("dirty agent triggers save on flush and updates state", {
  store    <- make_fake_store()
  tick     <- shiny::reactiveVal(0L)
  agent    <- make_stub_agent()
  agent_rv <- shiny::reactiveVal(agent)

  saved_agent <- make_saved_agent_for(agent)

  save_calls <- 0L
  pruned     <- FALSE

  local_mocked_bindings(
    session_is_dirty = function(session) TRUE,
    session_save = function(store, agent, ...) {
      save_calls <<- save_calls + 1L
      saved_agent
    },
    ovi_sessions = function(session_dir = NULL) {
      pruned <<- TRUE
      data.frame(session_id = character(0), updated_at = numeric(0),
                 stringsAsFactors = FALSE)
    },
    ovi_session_delete = function(session_id, session_dir = NULL) invisible(TRUE),
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      save_ctrl <<- copilot_save_controller(
        store                     = store,
        agent                     = agent_rv,
        history_invalidation_tick = tick,
        session                   = session
      )
    },
    expr = {
      save_ctrl$on_run_settled()
      # Before flush: scheduled, status flipped to "saving"
      expect_equal(shiny::isolate(save_ctrl$status()), "saving")
      expect_equal(save_calls, 0L)

      session$flushReact()

      expect_equal(save_calls, 1L)
      expect_true(pruned)
      expect_equal(shiny::isolate(save_ctrl$status()), "idle")
      expect_equal(shiny::isolate(tick()), 1L)
      # agent_rv was updated to the saved agent (same session_id)
      expect_equal(
        shiny::isolate(agent_rv())@session@saved_generation,
        saved_agent@session@saved_generation
      )
    }
  )
})

# ===========================================================================
# on_run_settled while save is already pending -> skip
# ===========================================================================

test_that("second on_run_settled while save is pending is skipped", {
  store    <- make_fake_store()
  tick     <- shiny::reactiveVal(0L)
  agent    <- make_stub_agent()
  agent_rv <- shiny::reactiveVal(agent)
  saved_agent <- make_saved_agent_for(agent)

  save_calls <- 0L
  local_mocked_bindings(
    session_is_dirty = function(session) TRUE,
    session_save = function(store, agent, ...) {
      save_calls <<- save_calls + 1L
      saved_agent
    },
    ovi_sessions = function(session_dir = NULL)
      data.frame(session_id = character(0), updated_at = numeric(0),
                 stringsAsFactors = FALSE),
    ovi_session_delete = function(session_id, session_dir = NULL) invisible(TRUE),
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      save_ctrl <<- copilot_save_controller(
        store                     = store,
        agent                     = agent_rv,
        history_invalidation_tick = tick,
        session                   = session
      )
    },
    expr = {
      save_ctrl$on_run_settled()
      expect_equal(shiny::isolate(save_ctrl$status()), "saving")
      # Second call while still "saving" must be a no-op
      expect_no_error(save_ctrl$on_run_settled())

      session$flushReact()
      # Only one save body executed
      expect_equal(save_calls, 1L)
      expect_equal(shiny::isolate(save_ctrl$status()), "idle")
    }
  )
})

# ===========================================================================
# session_id mismatch -> agent_rv NOT overwritten, but tick still bumped
# ===========================================================================

test_that("save result with mismatched session_id does not overwrite active agent", {
  store    <- make_fake_store()
  tick     <- shiny::reactiveVal(0L)
  original <- make_stub_agent()
  agent_rv <- shiny::reactiveVal(original)

  # session_save returns an agent with a different session_id (simulates
  # tier change between scheduling and onFlushed firing)
  diff_session <- S7::set_props(original@session, session_id = "different-id")
  diff_agent   <- S7::set_props(original, session = diff_session)
  original_sid <- original@session@session_id

  local_mocked_bindings(
    session_is_dirty = function(session) TRUE,
    session_save = function(store, agent, ...) diff_agent,
    ovi_sessions = function(session_dir = NULL)
      data.frame(session_id = character(0), updated_at = numeric(0),
                 stringsAsFactors = FALSE),
    ovi_session_delete = function(session_id, session_dir = NULL) invisible(TRUE),
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      save_ctrl <<- copilot_save_controller(
        store                     = store,
        agent                     = agent_rv,
        history_invalidation_tick = tick,
        session                   = session
      )
    },
    expr = {
      save_ctrl$on_run_settled()
      session$flushReact()
      # agent_rv keeps the original (session_id unchanged)
      expect_equal(shiny::isolate(agent_rv())@session@session_id, original_sid)
      # Tick still bumped — the saved session counts as persisted
      expect_equal(shiny::isolate(tick()), 1L)
      expect_equal(shiny::isolate(save_ctrl$status()), "idle")
    }
  )
})

# ===========================================================================
# session_save error -> status becomes "failed", no error propagates
# ===========================================================================

test_that("session_save failure flips status to 'failed' and is contained", {
  store    <- make_fake_store()
  tick     <- shiny::reactiveVal(0L)
  agent    <- make_stub_agent()
  agent_rv <- shiny::reactiveVal(agent)

  local_mocked_bindings(
    session_is_dirty = function(session) TRUE,
    session_save = function(store, agent, ...) stop("boom"),
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      save_ctrl <<- copilot_save_controller(
        store                     = store,
        agent                     = agent_rv,
        history_invalidation_tick = tick,
        session                   = session
      )
    },
    expr = {
      save_ctrl$on_run_settled()
      expect_no_error(session$flushReact())
      expect_equal(shiny::isolate(save_ctrl$status()), "failed")
      expect_equal(shiny::isolate(tick()), 0L)
    }
  )
})

# ===========================================================================
# Integration: real SessionStore on disk + flush plumbing
# ===========================================================================
#
# Confirms that the controller's onFlushed-deferred path drives a real
# SessionStore end-to-end. session_save is mocked (the package refuses to
# persist empty stub agents, and building a non-empty Agent in a unit test
# requires a live LLM round-trip we don't want here), but every other part
# of the path — onFlushed scheduling, status transitions, prune via real
# ovi_sessions on the real SQLite store, tick — is exercised against the
# real store and real Shiny session.

test_that("integration: real SessionStore + flush drives end-to-end save path", {
  tmpdir <- tempfile("copilot-save-int-")
  dir.create(tmpdir, recursive = TRUE)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  store <- omicsagentovi::SessionStore(session_dir = tmpdir)
  agent <- make_stub_agent()
  # Mark dirty so the save path actually runs
  dirty_session <- S7::set_props(
    agent@session,
    dirty_generation = agent@session@dirty_generation + 1L
  )
  agent <- S7::set_props(agent, session = dirty_session)
  expect_true(omicsagentovi::session_is_dirty(agent@session))

  saved_agent <- make_saved_agent_for(agent)

  local_mocked_bindings(
    session_save = function(store, agent, ...) saved_agent,
    .package = "omicsagentovi"
  )

  agent_rv <- shiny::reactiveVal(agent)
  tick     <- shiny::reactiveVal(0L)

  shiny::testServer(
    function(input, output, session) {
      save_ctrl <<- copilot_save_controller(
        store                     = store,
        agent                     = agent_rv,
        history_invalidation_tick = tick,
        session                   = session
      )
    },
    expr = {
      save_ctrl$on_run_settled()
      session$flushReact()

      after <- shiny::isolate(agent_rv())
      expect_false(omicsagentovi::session_is_dirty(after@session))
      expect_equal(shiny::isolate(save_ctrl$status()), "idle")
      expect_equal(shiny::isolate(tick()), 1L)
    }
  )
})
