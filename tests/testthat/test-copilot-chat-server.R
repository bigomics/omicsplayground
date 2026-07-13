## test-copilot-chat-server.R
##
## Tests for CopilotChatServer — chat module that owns shinychat I/O.
## Matches the test plan in .active_plans/refactor_copilot/chat/specs.md §"Test plan".

.ovi_source_dir <- Sys.getenv("OVI_SOURCE_DIR", "/Users/santiago/projects/BigOmics/omicsagentovi")
if (packageVersion("omicsagentovi") < "0.4.0" && dir.exists(.ovi_source_dir)) {
  pkgload::load_all(.ovi_source_dir, quiet = TRUE)
}

if (packageVersion("omicsagentovi") < "0.4.0") {
  skip("omicsagentovi >= 0.4.0 required.")
}

.board_dir <- if (dir.exists("components/app_copilot/R")) {
  "components/app_copilot/R"
} else {
  "../../components/app_copilot/R"
}

# Standalone bootstrap: %||% is normally in playbase / utils.R but is not on
# the search path when this file is sourced standalone in tests.
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

source(file.path(.board_dir, "copilot_options.R"),  local = TRUE)
source(file.path(.board_dir, "copilot_messages.R"), local = TRUE)
source(file.path(.board_dir, "copilot_logger.R"),   local = TRUE)
source(file.path(.board_dir, "CopilotChatServer.R"), local = TRUE)

library(omicsagentovi)

make_record <- function(role, content) {
  omicsagentovi::TranscriptRecord(
    idx = 1L, role = role, content_text = content, visibility = "visible"
  )
}

# ===========================================================================
# user input
# ===========================================================================

test_that("user text input triggers on_user_message with trimmed text", {
  captured <- NULL
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) captured <<- text,
      chat_event      = shiny::reactive(NULL)
    ), {
      session$setInputs(chat_user_input = "  hello world  ")
      expect_equal(captured, "hello world")
    }
  )
})

test_that("empty user input is ignored (no callback)", {
  captured <- NULL
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) captured <<- text,
      chat_event      = shiny::reactive(NULL)
    ), {
      session$setInputs(chat_user_input = "   ")
      expect_null(captured)
    }
  )
})

# ===========================================================================
# chat_event observer — post / clear / replay
# ===========================================================================

test_that("post event calls chat_append_message with the right role + content", {
  appended <- list()
  local_mocked_bindings(
    chat_append_message = function(id, message, chunk = FALSE) {
      appended[[length(appended) + 1L]] <<- message
      invisible(NULL)
    },
    chat_clear  = function(id) invisible(NULL),
    chat_append = function(id, ...) invisible(NULL),
    .package = "shinychat"
  )

  ev_rv <- shiny::reactiveVal(NULL)
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = ev_rv
    ), {
      ev_rv(list(type = "post", role = "assistant", text = "Hello there"))
      session$flushReact()

      expect_equal(length(appended), 1L)
      expect_equal(appended[[1]]$role, "assistant")
      expect_equal(appended[[1]]$content, "Hello there")
    }
  )
})

test_that("clear event calls chat_clear", {
  cleared <- 0L
  local_mocked_bindings(
    chat_clear          = function(id) { cleared <<- cleared + 1L; invisible(NULL) },
    chat_append_message = function(id, message, chunk = FALSE) invisible(NULL),
    chat_append         = function(id, ...) invisible(NULL),
    .package = "shinychat"
  )

  ev_rv <- shiny::reactiveVal(NULL)
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = ev_rv
    ), {
      ev_rv(list(type = "clear"))
      session$flushReact()
      expect_equal(cleared, 1L)
    }
  )
})

test_that("reset event calls chat_clear once then chat_append_message with greeting", {
  cleared  <- 0L
  appended <- list()
  local_mocked_bindings(
    chat_clear = function(id) { cleared <<- cleared + 1L; invisible(NULL) },
    chat_append_message = function(id, message, chunk = FALSE) {
      appended[[length(appended) + 1L]] <<- message
      invisible(NULL)
    },
    chat_append = function(id, ...) invisible(NULL),
    .package = "shinychat"
  )

  ev_rv <- shiny::reactiveVal(NULL)
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = ev_rv
    ), {
      ev_rv(list(type = "reset", role = "assistant", text = "Welcome!"))
      session$flushReact()

      expect_equal(cleared, 1L)
      expect_equal(length(appended), 1L)
      expect_equal(appended[[1]]$role,    "assistant")
      expect_equal(appended[[1]]$content, "Welcome!")
    }
  )
})

test_that("reset event with empty text only calls chat_clear (no post)", {
  cleared  <- 0L
  appended <- list()
  local_mocked_bindings(
    chat_clear = function(id) { cleared <<- cleared + 1L; invisible(NULL) },
    chat_append_message = function(id, message, chunk = FALSE) {
      appended[[length(appended) + 1L]] <<- message
      invisible(NULL)
    },
    chat_append = function(id, ...) invisible(NULL),
    .package = "shinychat"
  )

  ev_rv <- shiny::reactiveVal(NULL)
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = ev_rv
    ), {
      ev_rv(list(type = "reset", role = "assistant", text = ""))
      session$flushReact()

      expect_equal(cleared, 1L)
      expect_equal(length(appended), 0L)
    }
  )
})

test_that("replay event iterates records via chat_append_message", {
  appended <- list()
  local_mocked_bindings(
    chat_append_message = function(id, message, chunk = FALSE) {
      appended[[length(appended) + 1L]] <<- message
      invisible(NULL)
    },
    chat_clear  = function(id) invisible(NULL),
    chat_append = function(id, ...) invisible(NULL),
    .package = "shinychat"
  )

  recs <- list(
    make_record("user",      "Q1"),
    make_record("assistant", "A1"),
    make_record("user",      "Q2")
  )

  ev_rv <- shiny::reactiveVal(NULL)
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = ev_rv
    ), {
      ev_rv(list(type = "replay", records = recs))
      session$flushReact()

      expect_equal(length(appended), 3L)
      expect_equal(appended[[1]]$role,    "user")
      expect_equal(appended[[1]]$content, "Q1")
      expect_equal(appended[[2]]$role,    "assistant")
      expect_equal(appended[[3]]$content, "Q2")
    }
  )
})

# ===========================================================================
# stream event — promise resolution ticks stream_done
# ===========================================================================

test_that("stream event chains promise; stream_done ticks on resolve", {
  local_mocked_bindings(
    chat_append = function(id, gen, ...) promises::promise_resolve("ok"),
    chat_append_message = function(id, message, chunk = FALSE) invisible(NULL),
    chat_clear  = function(id) invisible(NULL),
    .package = "shinychat"
  )

  ev_rv <- shiny::reactiveVal(NULL)
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = ev_rv
    ), {
      api <- session$returned
      expect_equal(api$stream_done(), 0L)
      ev_rv(list(type = "stream", async_gen = list()))
      session$flushReact()
      # Drain the promises microtask queue.
      later::run_now()
      session$flushReact()
      expect_equal(api$stream_done(), 1L)
    }
  )
})

# ===========================================================================
# on_tool_request renders inline collapsible marker
# ===========================================================================

test_that("on_tool_request appends a collapsible <details> marker", {
  captured_content <- NULL
  local_mocked_bindings(
    chat_append_message = function(id, message, chunk = FALSE) {
      captured_content <<- message$content
      invisible(NULL)
    },
    chat_clear  = function(id) invisible(NULL),
    chat_append = function(id, ...) invisible(NULL),
    .package = "shinychat"
  )

  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = shiny::reactive(NULL)
    ), {
      api <- session$returned
      api$on_tool_request(list(
        name = "show_plot",
        arguments = list(plot_type = "pca"),
        id = "tool-123"
      ))
      expect_match(captured_content, "<details>")
      expect_match(captured_content, "show_plot")
      expect_match(captured_content, "pca")
    }
  )
})

# ===========================================================================
# tier_choices — updateRadioButtons is called with the right choices
# ===========================================================================

test_that("tier_choices update calls updateRadioButtons with correct choices", {
  calls <- list()
  local_mocked_bindings(
    updateRadioButtons = function(session, inputId, ...) {
      calls[[length(calls) + 1L]] <<- list(inputId = inputId, ...)
      invisible(NULL)
    },
    .package = "shiny"
  )

  choices_rv <- shiny::reactiveVal(NULL)
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = shiny::reactive(NULL),
      tier_choices    = choices_rv
    ), {
      ch <- c("Balanced (Default)" = "copilot-default", "Fast" = "copilot-fast")
      choices_rv(ch)
      session$flushReact()

      expect_true(length(calls) >= 1L)
      last_call <- calls[[length(calls)]]
      expect_equal(last_call$inputId, "tier_choice")
      expect_equal(last_call$choiceValues, unname(ch))
      expect_equal(last_call$choiceNames, names(ch))
    }
  )
})

# ===========================================================================
# tier_clicked — picking a tier via radioButtons updates the reactiveVal
# ===========================================================================

test_that("input$tier_choice triggers tier_clicked reactiveVal", {
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = shiny::reactive(NULL)
    ), {
      api <- session$returned
      expect_null(api$tier_clicked())

      # observeEvent with ignoreInit = TRUE skips the first setInputs call
      # (treated as init), so we need two distinct setInputs calls.
      session$setInputs(tier_choice = "copilot-default")
      session$flushReact()
      expect_null(api$tier_clicked())  # first call is the "init", ignored

      session$setInputs(tier_choice = "copilot-fast")
      session$flushReact()
      expect_equal(api$tier_clicked(), "copilot-fast")
    }
  )
})

# ===========================================================================
# starters — pushed inline as an assistant suggestion bubble after greeting
# ===========================================================================

test_that("first chat_event post triggers a starter suggestion bubble", {
  appended <- list()
  local_mocked_bindings(
    chat_append_message = function(id, message, chunk = FALSE) {
      appended[[length(appended) + 1L]] <<- message
      invisible(NULL)
    },
    chat_clear  = function(id) invisible(NULL),
    chat_append = function(id, ...) invisible(NULL),
    .package = "shinychat"
  )

  ev_rv <- shiny::reactiveVal(NULL)
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = ev_rv,
      starters        = shiny::reactive(c("Summarize findings", "Top biomarkers"))
    ), {
      ev_rv(list(type = "post", role = "assistant", text = "Hi"))
      session$flushReact()
      # Expect 2 messages: greeting + starter bubble.
      expect_equal(length(appended), 2L)
      expect_equal(appended[[1]]$content, "Hi")
      expect_match(appended[[2]]$content, "copilot-starter-list", fixed = TRUE)
      expect_match(appended[[2]]$content, "suggestion submit",    fixed = TRUE)
      expect_match(appended[[2]]$content, "Summarize findings",   fixed = TRUE)
      expect_match(appended[[2]]$content, "Top biomarkers",       fixed = TRUE)
    }
  )
})

test_that("reset event re-seeds the starter bubble", {
  appended <- list()
  local_mocked_bindings(
    chat_append_message = function(id, message, chunk = FALSE) {
      appended[[length(appended) + 1L]] <<- message
      invisible(NULL)
    },
    chat_clear  = function(id) invisible(NULL),
    chat_append = function(id, ...) invisible(NULL),
    .package = "shinychat"
  )

  ev_rv <- shiny::reactiveVal(NULL)
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = ev_rv,
      starters        = shiny::reactive(c("Q1", "Q2"))
    ), {
      ev_rv(list(type = "reset", role = "assistant", text = "Welcome"))
      session$flushReact()
      # Expect 2 messages: greeting + starter bubble.
      expect_equal(length(appended), 2L)
      expect_equal(appended[[1]]$content, "Welcome")
      expect_match(appended[[2]]$content, "Q1", fixed = TRUE)
      expect_match(appended[[2]]$content, "Q2", fixed = TRUE)
    }
  )
})

test_that("subsequent posts after the first do not re-append starters", {
  appended <- list()
  local_mocked_bindings(
    chat_append_message = function(id, message, chunk = FALSE) {
      appended[[length(appended) + 1L]] <<- message
      invisible(NULL)
    },
    chat_clear  = function(id) invisible(NULL),
    chat_append = function(id, ...) invisible(NULL),
    .package = "shinychat"
  )

  ev_rv <- shiny::reactiveVal(NULL)
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = ev_rv,
      starters        = shiny::reactive(c("Q1", "Q2"))
    ), {
      ev_rv(list(type = "post", role = "assistant", text = "Hi"))
      session$flushReact()
      ev_rv(list(type = "post", role = "user", text = "hello"))
      session$flushReact()
      ev_rv(list(type = "post", role = "assistant", text = "reply"))
      session$flushReact()
      # Greeting + starter bubble + 2 follow-up posts = 4 messages.
      expect_equal(length(appended), 4L)
      # No additional starter bubble on subsequent posts.
      expect_false(any(grepl("copilot-starter-list",
                              vapply(appended[3:4], function(m) m$content, character(1)),
                              fixed = TRUE)))
    }
  )
})

# ===========================================================================
# Public surface shape
# ===========================================================================

test_that("CopilotChatServer returns the documented API list", {
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = shiny::reactive(NULL)
    ), {
      api <- session$returned
      expect_true(shiny::is.reactive(api$user_input))
      expect_true(shiny::is.reactive(api$stream_done))
      expect_true(shiny::is.reactive(api$last_error))
      expect_true(is.function(api$on_tool_request))
      expect_null(api$push_event)
      # tier_clicked is a reactiveVal (also a function/reactive)
      expect_true(is.function(api$tier_clicked))
    }
  )
})

# ===========================================================================
# Bug 1 — stream rejection calls on_stream_error and ticks stream_done
# ===========================================================================

test_that("stream rejection calls on_stream_error and ticks stream_done", {
  stream_error_called <- FALSE
  appended <- list()
  local_mocked_bindings(
    chat_append = function(id, gen, ...) promises::promise_reject(
      simpleError("insufficient_quota: you exceeded your current quota")
    ),
    chat_append_message = function(id, message, chunk = FALSE) {
      appended[[length(appended) + 1L]] <<- message
      invisible(NULL)
    },
    chat_clear  = function(id) invisible(NULL),
    .package = "shinychat"
  )

  ev_rv <- shiny::reactiveVal(NULL)
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = ev_rv,
      on_stream_error = function() { stream_error_called <<- TRUE }
    ), {
      api <- session$returned
      expect_equal(api$stream_done(), 0L)

      ev_rv(list(type = "stream", async_gen = list()))
      session$flushReact()
      later::run_now()
      session$flushReact()

      # on_stream_error must have been called to settle the run lifecycle
      expect_true(stream_error_called)
      # stream_done ticks so run_controller's observer can fire
      expect_equal(api$stream_done(), 1L)
      # last_error captures the raw message for debugging
      expect_match(api$last_error(), "quota", fixed = TRUE)
    }
  )
})

test_that("stream rejection with no on_stream_error still ticks stream_done", {
  local_mocked_bindings(
    chat_append = function(id, gen, ...) promises::promise_reject(
      simpleError("some error")
    ),
    chat_append_message = function(id, message, chunk = FALSE) invisible(NULL),
    chat_clear  = function(id) invisible(NULL),
    .package = "shinychat"
  )

  ev_rv <- shiny::reactiveVal(NULL)
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = ev_rv
      # on_stream_error omitted — defaults to NULL
    ), {
      api <- session$returned
      ev_rv(list(type = "stream", async_gen = list()))
      session$flushReact()
      later::run_now()
      session$flushReact()
      expect_equal(api$stream_done(), 1L)
    }
  )
})

# ===========================================================================
# Bug 2 — stream rejection posts a friendly message, not raw API text
# ===========================================================================

test_that("stream rejection posts a friendly quota message for known quota error", {
  posted_content <- NULL
  local_mocked_bindings(
    chat_append = function(id, gen, ...) promises::promise_reject(
      simpleError("openai: insufficient_quota — you exceeded your monthly quota")
    ),
    chat_append_message = function(id, message, chunk = FALSE) {
      posted_content <<- message$content
      invisible(NULL)
    },
    chat_clear  = function(id) invisible(NULL),
    .package = "shinychat"
  )

  ev_rv <- shiny::reactiveVal(NULL)
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = ev_rv
    ), {
      ev_rv(list(type = "stream", async_gen = list()))
      session$flushReact()
      later::run_now()
      session$flushReact()

      # Must NOT contain raw API error text
      expect_false(grepl("insufficient_quota", posted_content, fixed = TRUE))
      # Must contain billing/quota guidance
      expect_match(posted_content, "quota", fixed = TRUE)
    }
  )
})

test_that("stream rejection posts a friendly auth message for 401/invalid key", {
  posted_content <- NULL
  local_mocked_bindings(
    chat_append = function(id, gen, ...) promises::promise_reject(
      simpleError("401 Unauthorized: invalid api key provided")
    ),
    chat_append_message = function(id, message, chunk = FALSE) {
      posted_content <<- message$content
      invisible(NULL)
    },
    chat_clear  = function(id) invisible(NULL),
    .package = "shinychat"
  )

  ev_rv <- shiny::reactiveVal(NULL)
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = ev_rv
    ), {
      ev_rv(list(type = "stream", async_gen = list()))
      session$flushReact()
      later::run_now()
      session$flushReact()

      expect_false(grepl("invalid api key", posted_content, fixed = TRUE))
      expect_match(posted_content, "API key", fixed = TRUE)
    }
  )
})

test_that("stream rejection posts a generic message for unrecognised errors", {
  posted_content <- NULL
  local_mocked_bindings(
    chat_append = function(id, gen, ...) promises::promise_reject(
      simpleError("something weird happened with the model backend xyz-123")
    ),
    chat_append_message = function(id, message, chunk = FALSE) {
      posted_content <<- message$content
      invisible(NULL)
    },
    chat_clear  = function(id) invisible(NULL),
    .package = "shinychat"
  )

  ev_rv <- shiny::reactiveVal(NULL)
  shiny::testServer(CopilotChatServer,
    args = list(
      id              = "chat",
      on_user_message = function(text) NULL,
      chat_event      = ev_rv
    ), {
      ev_rv(list(type = "stream", async_gen = list()))
      session$flushReact()
      later::run_now()
      session$flushReact()

      expect_false(grepl("xyz-123", posted_content, fixed = TRUE))
      expect_match(posted_content, "couldn't complete", fixed = TRUE)
    }
  )
})

# ===========================================================================
# copilot_friendly_error — message mapping unit tests
# ===========================================================================

test_that("copilot_friendly_error maps quota errors correctly", {
  expect_match(copilot_friendly_error("insufficient_quota: monthly cap reached"), "quota")
  expect_match(copilot_friendly_error("you exceeded your current quota"), "quota")
  expect_match(copilot_friendly_error("429 Too Many Requests: billing limit"), "quota")
})

test_that("copilot_friendly_error maps auth errors correctly", {
  expect_match(copilot_friendly_error("401 Unauthorized"), "API key")
  expect_match(copilot_friendly_error("invalid api key provided"), "API key")
  expect_match(copilot_friendly_error("authentication failed"), "API key")
})

test_that("copilot_friendly_error maps rate-limit errors correctly", {
  expect_match(copilot_friendly_error("429 Too Many Requests: rate limit exceeded"), "rate")
  expect_match(copilot_friendly_error("rate_limit: too many requests"), "rate")
})

test_that("copilot_friendly_error returns generic message for unknown errors", {
  result <- copilot_friendly_error("an internal model error occurred")
  expect_match(result, "couldn't complete", fixed = TRUE)
})

test_that("copilot_friendly_error accepts a condition object", {
  e <- simpleError("insufficient_quota: over limit")
  expect_match(copilot_friendly_error(e), "quota")
})
