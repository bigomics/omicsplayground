##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# Bootstrap: locate component root from test runner or IDE
.module_dir <- if (dir.exists("components/modules")) {
  "components/modules"
} else {
  "../../components/modules"
}

# Minimal stub for getUserOption — mirrors utils.R:79
getUserOption <- function(session, var, value) session$userData[[var]]

source(file.path(.module_dir, "AiCards.R"), local = TRUE)

# A lightweight fake session that carries userData like the real Shiny session.
make_fake_session <- function(userData = list()) {
  list(userData = list2env(userData, parent = emptyenv()))
}

test_that("get_ai_credentials returns NULL when nothing is set", {
  session <- make_fake_session()
  expect_null(get_ai_credentials(session))
})

test_that("get_ai_credentials returns NULL for bigomics provider (key stored as NULL)", {
  session <- make_fake_session(list(ai_credentials = NULL))
  expect_null(get_ai_credentials(session))
})

test_that("get_ai_credentials returns the closure when a BYOK key is stored", {
  key <- "sk-test-1234"
  cred_fn <- local({ k <- key; function() k })
  session <- make_fake_session(list(ai_credentials = cred_fn))

  result <- get_ai_credentials(session)
  expect_true(is.function(result))
  expect_equal(result(), "sk-test-1234")
})

test_that("get_ai_credentials closure captures the value, not the binding", {
  key <- "original-key"
  cred_fn <- local({ k <- key; function() k })
  session <- make_fake_session(list(ai_credentials = cred_fn))

  # Even if the outer variable changes, the captured value stays fixed.
  key <- "changed-key"
  result <- get_ai_credentials(session)
  expect_equal(result(), "original-key")
})

test_that("get_ai_credentials result is suitable for future_promise capture", {
  # Verify the closure is a plain R object (not a reactive), safe to
  # serialize into a future worker.
  key <- "byok-key"
  cred_fn <- local({ k <- key; function() k })
  session <- make_fake_session(list(ai_credentials = cred_fn))

  captured <- get_ai_credentials(session)
  # A plain function — not a reactive, not a promise.
  expect_true(is.function(captured))
  expect_false(inherits(captured, "reactive"))
  expect_equal(captured(), "byok-key")
})
