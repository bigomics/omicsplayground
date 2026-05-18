## test-copilot-run-controller.R
##
## Tests for copilot_run_controller() — unified dispatch for ask/new_chat/
## tier/abort. Matches the test plan in
## .active_plans/refactor_copilot/run_controller/specs.md §"Test plan".

.ovi_source_dir <- Sys.getenv("OVI_SOURCE_DIR", "/Users/santiago/projects/BigOmics/omicsagentovi")
if (packageVersion("omicsagentovi") < "0.4.0" && dir.exists(.ovi_source_dir)) {
  pkgload::load_all(.ovi_source_dir, quiet = TRUE)
}

if (packageVersion("omicsagentovi") < "0.4.0") {
  skip("omicsagentovi >= 0.4.0 required.")
}

.board_dir <- if (dir.exists("components/board.copilot/R")) {
  "components/board.copilot/R"
} else {
  "../../components/board.copilot/R"
}

source(file.path(.board_dir, "copilot_options.R"),  local = TRUE)
source(file.path(.board_dir, "copilot_messages.R"), local = TRUE)
source(file.path(.board_dir, "copilot_logger.R"),   local = TRUE)
source(file.path(.board_dir, "copilot_bindings.R"),          local = TRUE)
source(file.path(.board_dir, "copilot_run_controller.R"),    local = TRUE)

library(omicsagentovi)

.make_stub_chat <- function() {
  list(get_turns = function() list(), set_turns = function(x) invisible(NULL))
}
make_stub_agent <- function() {
  omicsagentovi::Agent(tools = list(), chat = .make_stub_chat())
}
make_fake_store <- function() {
  omicsagentovi::SessionStore(session_dir = tempdir())
}

# ===========================================================================
# Request constructors
# ===========================================================================

test_that("run_request_ask returns kind=ask with text + show_user_msg", {
  r <- run_request_ask("hello")
  expect_equal(r$kind, "ask")
  expect_equal(r$text, "hello")
  expect_true(r$show_user_msg)

  r2 <- run_request_ask("hi", show_user_msg = FALSE)
  expect_false(r2$show_user_msg)
})

test_that("run_request_new_chat returns kind=new_chat", {
  expect_equal(run_request_new_chat()$kind, "new_chat")
})

test_that("run_request_tier returns kind=tier with new_tier", {
  r <- run_request_tier("copilot-deep")
  expect_equal(r$kind, "tier")
  expect_equal(r$new_tier, "copilot-deep")
})

test_that("run_request_abort returns kind=abort with reason", {
  r <- run_request_abort()
  expect_equal(r$kind, "abort")
  expect_match(r$reason, "stopped", fixed = TRUE)
})

# ===========================================================================
# Helper to construct a controller wired with mocks
# ===========================================================================

# Use within testServer: it provides `session`. The factory keeps closures
# over the reactives passed in.

# ===========================================================================
# dispatch guards
# ===========================================================================

test_that("dispatch returns early while streaming for non-abort kinds", {
  agent_rv      <- shiny::reactiveVal(make_stub_agent())
  run_status_rv <- shiny::reactiveVal("streaming")
  tier_rv       <- shiny::reactiveVal("copilot-default")
  pgx_evt       <- shiny::reactiveVal(NULL)
  chat_evt      <- shiny::reactiveVal(NULL)
  pgx           <- shiny::reactiveValues(name = NULL)

  prompt_stream_called <- FALSE
  local_mocked_bindings(
    agent_prompt_stream = function(...) {
      prompt_stream_called <<- TRUE
      list()
    },
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_run_controller(
        agent = agent_rv, run_status = run_status_rv, tier = tier_rv,
        save_ctrl = list(on_run_settled = function() NULL),
        restore_ctrl = NULL, evidence = NULL,
        store = make_fake_store(),
        chat_event = chat_evt,
        chat_on_tool_request = function(req) NULL,
        pgx = pgx, pgx_dir = tempdir(), docs_dir = tempdir(),
        pgx_loaded_event = pgx_evt, maxturns = Inf, session = session
      )
    },
    expr = {
      ctrl$dispatch(run_request_ask("hi"))
      expect_false(prompt_stream_called)
      expect_equal(shiny::isolate(agent_rv())@session@session_id,
                   shiny::isolate(agent_rv())@session@session_id)
    }
  )
})

test_that("dispatch(ask) with NULL agent posts a hint and does not stream", {
  agent_rv      <- shiny::reactiveVal(NULL)
  run_status_rv <- shiny::reactiveVal("idle")
  tier_rv       <- shiny::reactiveVal("copilot-default")
  pgx_evt       <- shiny::reactiveVal(NULL)
  chat_evt      <- shiny::reactiveVal(NULL)
  pgx           <- shiny::reactiveValues(name = NULL)

  posted <- list()
  push_writer <- function(ev) {
    if (missing(ev)) return(chat_evt())
    posted[[length(posted) + 1L]] <<- ev
    chat_evt(ev)
  }

  prompt_called <- FALSE
  local_mocked_bindings(
    agent_prompt_stream = function(...) { prompt_called <<- TRUE; list() },
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_run_controller(
        agent = agent_rv, run_status = run_status_rv, tier = tier_rv,
        save_ctrl = list(on_run_settled = function() NULL),
        restore_ctrl = NULL, evidence = NULL,
        store = make_fake_store(),
        chat_event = push_writer,
        chat_on_tool_request = function(req) NULL,
        pgx = pgx, pgx_dir = tempdir(), docs_dir = tempdir(),
        pgx_loaded_event = pgx_evt, maxturns = Inf, session = session
      )
    },
    expr = {
      ctrl$dispatch(run_request_ask("hi"))
      expect_false(prompt_called)
      expect_true(length(posted) >= 1L)
      expect_equal(posted[[1]]$type, "post")
      expect_match(posted[[1]]$text, "dataset", fixed = TRUE)
    }
  )
})

# ===========================================================================
# dispatch(ask) happy path
# ===========================================================================

test_that("dispatch(ask) happy path: user msg posted, stream invoked, status set", {
  agent_rv      <- shiny::reactiveVal(make_stub_agent())
  run_status_rv <- shiny::reactiveVal("idle")
  tier_rv       <- shiny::reactiveVal("copilot-default")
  pgx_evt       <- shiny::reactiveVal(NULL)
  chat_evt      <- shiny::reactiveVal(NULL)
  pgx           <- shiny::reactiveValues(name = NULL)

  posted <- list()
  push_writer <- function(ev) {
    if (missing(ev)) return(chat_evt())
    posted[[length(posted) + 1L]] <<- ev
    chat_evt(ev)
  }

  captured_callbacks <- list()
  sentinel <- structure(list(), class = "fake_async_gen")

  local_mocked_bindings(
    agent_prompt_stream = function(agent, text, ..., on_done = NULL,
                                   on_tool_request = NULL) {
      captured_callbacks$on_done           <<- on_done
      captured_callbacks$on_tool_request   <<- on_tool_request
      captured_callbacks$text              <<- text
      sentinel
    },
    session_is_dirty = function(s) FALSE,
    session_transcript = function(s, view = "user") list(),
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_run_controller(
        agent = agent_rv, run_status = run_status_rv, tier = tier_rv,
        save_ctrl = list(on_run_settled = function() NULL),
        restore_ctrl = NULL, evidence = NULL,
        store = make_fake_store(),
        chat_event = push_writer,
        chat_on_tool_request = function(req) "tool-cb",
        pgx = pgx, pgx_dir = tempdir(), docs_dir = tempdir(),
        pgx_loaded_event = pgx_evt, maxturns = Inf, session = session
      )
    },
    expr = {
      ctrl$dispatch(run_request_ask("what is up?", show_user_msg = TRUE))

      # Events: user post + stream
      expect_true(length(posted) >= 2L)
      kinds <- vapply(posted, function(e) e$type, character(1))
      expect_true("post"   %in% kinds)
      expect_true("stream" %in% kinds)

      # Run status flipped to streaming
      expect_equal(shiny::isolate(run_status_rv()), "streaming")

      # agent_prompt_stream received our callbacks
      expect_true(is.function(captured_callbacks$on_done))
      expect_true(is.function(captured_callbacks$on_tool_request))
      expect_equal(captured_callbacks$text, "what is up?")
    }
  )
})

test_that("on_done callback updates agent and run_status; calls save when dirty", {
  agent_rv      <- shiny::reactiveVal(make_stub_agent())
  run_status_rv <- shiny::reactiveVal("idle")
  tier_rv       <- shiny::reactiveVal("copilot-default")
  pgx_evt       <- shiny::reactiveVal(NULL)
  chat_evt      <- shiny::reactiveVal(NULL)
  pgx           <- shiny::reactiveValues(name = NULL)

  save_called <- FALSE
  save_ctrl <- list(on_run_settled = function() save_called <<- TRUE)

  captured_on_done <- NULL
  refreshed_agent  <- make_stub_agent()

  local_mocked_bindings(
    agent_prompt_stream = function(agent, text, ..., on_done = NULL,
                                   on_tool_request = NULL) {
      captured_on_done <<- on_done
      list()
    },
    session_is_dirty   = function(s) TRUE,
    session_transcript = function(s, view = "user") list(),
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_run_controller(
        agent = agent_rv, run_status = run_status_rv, tier = tier_rv,
        save_ctrl = save_ctrl, restore_ctrl = NULL, evidence = NULL,
        store = make_fake_store(),
        chat_event = chat_evt, chat_on_tool_request = function(req) NULL,
        pgx = pgx, pgx_dir = tempdir(), docs_dir = tempdir(),
        pgx_loaded_event = pgx_evt, maxturns = Inf, session = session
      )
    },
    expr = {
      ctrl$dispatch(run_request_ask("hi"))
      expect_true(is.function(captured_on_done))

      # Simulate stream settled
      captured_on_done(list(agent = refreshed_agent, status = "completed",
                            text = "x", error = NULL))
      expect_identical(shiny::isolate(agent_rv()), refreshed_agent)
      expect_equal(shiny::isolate(run_status_rv()), "completed")
      expect_true(save_called)
    }
  )
})

# ===========================================================================
# Abort
# ===========================================================================

test_that("dispatch(abort) calls agent_request_abort with reason", {
  agent_rv      <- shiny::reactiveVal(make_stub_agent())
  run_status_rv <- shiny::reactiveVal("streaming")
  tier_rv       <- shiny::reactiveVal("copilot-default")
  pgx_evt       <- shiny::reactiveVal(NULL)
  chat_evt      <- shiny::reactiveVal(NULL)
  pgx           <- shiny::reactiveValues(name = NULL)

  captured_reason <- NULL
  local_mocked_bindings(
    agent_request_abort = function(agent, reason = "User stopped the run") {
      captured_reason <<- reason
      TRUE
    },
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_run_controller(
        agent = agent_rv, run_status = run_status_rv, tier = tier_rv,
        save_ctrl = list(on_run_settled = function() NULL),
        restore_ctrl = NULL, evidence = NULL,
        store = make_fake_store(),
        chat_event = chat_evt, chat_on_tool_request = function(req) NULL,
        pgx = pgx, pgx_dir = tempdir(), docs_dir = tempdir(),
        pgx_loaded_event = pgx_evt, maxturns = Inf, session = session
      )
    },
    expr = {
      ctrl$dispatch(run_request_abort("stopit"))
      expect_equal(captured_reason, "stopit")
    }
  )
})

# ===========================================================================
# new_chat / tier — reset path
# ===========================================================================

test_that("dispatch(new_chat) with dirty agent pre-saves and rebuilds; clears chat", {
  agent_rv      <- shiny::reactiveVal(make_stub_agent())
  run_status_rv <- shiny::reactiveVal("idle")
  tier_rv       <- shiny::reactiveVal("copilot-default")
  pgx_evt       <- shiny::reactiveVal(NULL)
  chat_evt      <- shiny::reactiveVal(NULL)
  pgx           <- shiny::reactiveValues(name = NULL)

  posted <- list()
  push_writer <- function(ev) {
    if (missing(ev)) return(chat_evt())
    posted[[length(posted) + 1L]] <<- ev
    chat_evt(ev)
  }

  session_save_called <- FALSE
  saved_agent <- make_stub_agent()

  local_mocked_bindings(
    session_is_dirty = function(s) TRUE,
    session_save = function(store, agent) {
      session_save_called <<- TRUE
      saved_agent
    },
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_run_controller(
        agent = agent_rv, run_status = run_status_rv, tier = tier_rv,
        save_ctrl = list(on_run_settled = function() NULL),
        restore_ctrl = NULL, evidence = NULL,
        store = make_fake_store(),
        chat_event = push_writer, chat_on_tool_request = function(req) NULL,
        pgx = pgx, pgx_dir = tempdir(), docs_dir = tempdir(),
        pgx_loaded_event = pgx_evt, maxturns = Inf, session = session
      )
    },
    expr = {
      ctrl$dispatch(run_request_new_chat())

      expect_true(session_save_called)

      # With NULL pgx, the reset leaves agent as NULL (no Agent constructed)
      expect_null(shiny::isolate(agent_rv()))

      # Chat events: single reset event (clear + greeting in one flush-safe write)
      kinds <- vapply(posted, function(e) e$type, character(1))
      expect_true("reset" %in% kinds)

      expect_equal(shiny::isolate(run_status_rv()), "idle")
    }
  )
})

test_that("dispatch(tier) with same tier is a no-op", {
  agent_rv      <- shiny::reactiveVal(make_stub_agent())
  run_status_rv <- shiny::reactiveVal("idle")
  tier_rv       <- shiny::reactiveVal("copilot-default")
  pgx_evt       <- shiny::reactiveVal(NULL)
  chat_evt      <- shiny::reactiveVal(NULL)
  pgx           <- shiny::reactiveValues(name = NULL)

  session_save_called <- FALSE
  local_mocked_bindings(
    session_is_dirty = function(s) TRUE,
    session_save = function(store, agent) {
      session_save_called <<- TRUE
      agent
    },
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_run_controller(
        agent = agent_rv, run_status = run_status_rv, tier = tier_rv,
        save_ctrl = list(on_run_settled = function() NULL),
        restore_ctrl = NULL, evidence = NULL,
        store = make_fake_store(),
        chat_event = chat_evt, chat_on_tool_request = function(req) NULL,
        pgx = pgx, pgx_dir = tempdir(), docs_dir = tempdir(),
        pgx_loaded_event = pgx_evt, maxturns = Inf, session = session
      )
    },
    expr = {
      ctrl$dispatch(run_request_tier("copilot-default"))
      expect_false(session_save_called)
    }
  )
})

test_that("dispatch(tier) with different tier pre-saves dirty + updates tier", {
  agent_rv      <- shiny::reactiveVal(make_stub_agent())
  run_status_rv <- shiny::reactiveVal("idle")
  tier_rv       <- shiny::reactiveVal("copilot-default")
  pgx_evt       <- shiny::reactiveVal(NULL)
  chat_evt      <- shiny::reactiveVal(NULL)
  pgx           <- shiny::reactiveValues(name = NULL)

  session_save_called <- FALSE
  local_mocked_bindings(
    session_is_dirty = function(s) TRUE,
    session_save = function(store, agent) {
      session_save_called <<- TRUE
      agent
    },
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_run_controller(
        agent = agent_rv, run_status = run_status_rv, tier = tier_rv,
        save_ctrl = list(on_run_settled = function() NULL),
        restore_ctrl = NULL, evidence = NULL,
        store = make_fake_store(),
        chat_event = chat_evt, chat_on_tool_request = function(req) NULL,
        pgx = pgx, pgx_dir = tempdir(), docs_dir = tempdir(),
        pgx_loaded_event = pgx_evt, maxturns = Inf, session = session
      )
    },
    expr = {
      ctrl$dispatch(run_request_tier("copilot-deep"))
      expect_true(session_save_called)
      expect_equal(shiny::isolate(tier_rv()), "copilot-deep")
    }
  )
})

# ===========================================================================
# Maxturns guard
# ===========================================================================

test_that("dispatch(ask) refuses when transcript user-turns >= maxturns", {
  agent_rv      <- shiny::reactiveVal(make_stub_agent())
  run_status_rv <- shiny::reactiveVal("idle")
  tier_rv       <- shiny::reactiveVal("copilot-default")
  pgx_evt       <- shiny::reactiveVal(NULL)
  chat_evt      <- shiny::reactiveVal(NULL)
  pgx           <- shiny::reactiveValues(name = NULL)

  prompt_called <- FALSE
  fake_records <- list(
    omicsagentovi::TranscriptRecord(idx=1L, role="user",
      content_text="x", visibility="visible"),
    omicsagentovi::TranscriptRecord(idx=2L, role="user",
      content_text="y", visibility="visible")
  )

  local_mocked_bindings(
    agent_prompt_stream = function(...) { prompt_called <<- TRUE; list() },
    session_transcript  = function(s, view = "user") fake_records,
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_run_controller(
        agent = agent_rv, run_status = run_status_rv, tier = tier_rv,
        save_ctrl = list(on_run_settled = function() NULL),
        restore_ctrl = NULL, evidence = NULL,
        store = make_fake_store(),
        chat_event = chat_evt, chat_on_tool_request = function(req) NULL,
        pgx = pgx, pgx_dir = tempdir(), docs_dir = tempdir(),
        pgx_loaded_event = pgx_evt, maxturns = 2L, session = session
      )
    },
    expr = {
      ctrl$dispatch(run_request_ask("third question"))
      expect_false(prompt_called)
    }
  )
})

# ===========================================================================
# apply_dataset
# ===========================================================================

test_that("apply_dataset with NULL agent constructs a new Agent", {
  agent_rv      <- shiny::reactiveVal(NULL)
  run_status_rv <- shiny::reactiveVal("idle")
  tier_rv       <- shiny::reactiveVal("copilot-default")
  pgx_evt       <- shiny::reactiveVal(NULL)
  chat_evt      <- shiny::reactiveVal(NULL)
  pgx           <- shiny::reactiveValues(name = NULL)

  agent_constructed <- FALSE
  fake_agent <- make_stub_agent()

  local_mocked_bindings(
    Agent = function(tier, context, session, bindings) {
      agent_constructed <<- TRUE
      fake_agent
    },
    RunContext = function(pgx) list(pgx = pgx),
    AgentSession = function(session_id = NA_character_) list(session_id = session_id),
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_run_controller(
        agent = agent_rv, run_status = run_status_rv, tier = tier_rv,
        save_ctrl = list(on_run_settled = function() NULL),
        restore_ctrl = NULL, evidence = NULL,
        store = make_fake_store(),
        chat_event = chat_evt, chat_on_tool_request = function(req) NULL,
        pgx = pgx, pgx_dir = tempdir(), docs_dir = tempdir(),
        pgx_loaded_event = pgx_evt, maxturns = Inf, session = session
      )
    },
    expr = {
      ctrl$apply_dataset(list(name = "ds1"), "ds1", NULL, tempdir())
      expect_true(agent_constructed)
      expect_identical(shiny::isolate(agent_rv()), fake_agent)
    }
  )
})

test_that("apply_dataset with existing agent calls agent_set_pgx (chat NOT cleared)", {
  initial_agent <- make_stub_agent()
  updated_agent <- make_stub_agent()
  agent_rv      <- shiny::reactiveVal(initial_agent)
  run_status_rv <- shiny::reactiveVal("idle")
  tier_rv       <- shiny::reactiveVal("copilot-default")
  pgx_evt       <- shiny::reactiveVal(NULL)
  chat_evt      <- shiny::reactiveVal(NULL)
  pgx           <- shiny::reactiveValues(name = NULL)

  posted <- list()
  push_writer <- function(ev) {
    if (missing(ev)) return(chat_evt())
    posted[[length(posted) + 1L]] <<- ev
    chat_evt(ev)
  }

  set_pgx_called <- FALSE
  local_mocked_bindings(
    agent_set_pgx = function(agent, pgx, ...) {
      set_pgx_called <<- TRUE
      updated_agent
    },
    .package = "omicsagentovi"
  )

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_run_controller(
        agent = agent_rv, run_status = run_status_rv, tier = tier_rv,
        save_ctrl = list(on_run_settled = function() NULL),
        restore_ctrl = NULL, evidence = NULL,
        store = make_fake_store(),
        chat_event = push_writer, chat_on_tool_request = function(req) NULL,
        pgx = pgx, pgx_dir = tempdir(), docs_dir = tempdir(),
        pgx_loaded_event = pgx_evt, maxturns = Inf, session = session
      )
    },
    expr = {
      ctrl$apply_dataset(list(name = "ds2"), "ds2", NULL, tempdir())
      expect_true(set_pgx_called)
      expect_identical(shiny::isolate(agent_rv()), updated_agent)

      # Chat preserved — no clear event pushed by apply_dataset
      kinds <- vapply(posted, function(e) e$type, character(1))
      expect_false("clear" %in% kinds)
    }
  )
})

# ===========================================================================
# Public surface shape
# ===========================================================================

test_that("copilot_run_controller returns dispatch + apply_dataset", {
  agent_rv      <- shiny::reactiveVal(NULL)
  run_status_rv <- shiny::reactiveVal("idle")
  tier_rv       <- shiny::reactiveVal("copilot-default")
  pgx_evt       <- shiny::reactiveVal(NULL)
  chat_evt      <- shiny::reactiveVal(NULL)
  pgx           <- shiny::reactiveValues(name = NULL)

  shiny::testServer(
    function(input, output, session) {
      ctrl <<- copilot_run_controller(
        agent = agent_rv, run_status = run_status_rv, tier = tier_rv,
        save_ctrl = list(on_run_settled = function() NULL),
        restore_ctrl = NULL, evidence = NULL,
        store = make_fake_store(),
        chat_event = chat_evt, chat_on_tool_request = function(req) NULL,
        pgx = pgx, pgx_dir = tempdir(), docs_dir = tempdir(),
        pgx_loaded_event = pgx_evt, maxturns = Inf, session = session
      )
    },
    expr = {
      expect_true(is.function(ctrl$dispatch))
      expect_true(is.function(ctrl$apply_dataset))
    }
  )
})
