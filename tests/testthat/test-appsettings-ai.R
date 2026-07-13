## test-appsettings-ai.R
##
## Unit tests for the AI Features server logic in AppSettingsBoard
## (Task 2.3 of the BYOK AI-providers feature). Exercises the provider /
## credential store, provider-aware menu repopulation, the admin lock and
## the logout clear via shiny::testServer + testthat mocked bindings.

suppressMessages(library(shiny))

## ---- test bootstrap: make the module's free variables resolvable ----------

## Standard null-coalescing operator (test bootstrap only; never redefine in R/).
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

## session$userData-backed option store (mirrors app/R/utils/utils.R).
setUserOption <- function(session, var, value) session$userData[[var]] <- value
getUserOption <- function(session, var, value) session$userData[[var]]

## Init-time collaborators the module touches but that are irrelevant here.
dbg                         <- function(...) invisible(NULL)
user_table_resources_server <- function(...) invisible(NULL)
get_color_theme             <- function() shiny::reactiveValues()
load_color_theme            <- function(...) NULL
save_color_theme            <- function(...) invisible(NULL)
PlotModuleServer            <- function(...) invisible(NULL)
TableModuleServer           <- function(...) invisible(NULL)
COLOR_THEME_DEFAULTS        <- list()

.repo_dir <- if (dir.exists("components/app/R")) {
  normalizePath(".")
} else {
  normalizePath("../..")
}
.required_omicsai_exports <- c(
  "ai_known_models",
  "ai_provider_catalog",
  "ai_select_model",
  "ai_validate_model",
  "ai.list_provider_models"
)
testthat::skip_if_not_installed("omicsai", minimum_version = "0.3.2")
testthat::skip_if_not(
  all(.required_omicsai_exports %in% getNamespaceExports("omicsai")),
  "omicsai provider catalog API exports are required"
)
source(file.path(.repo_dir, "components/app/R/ai_model_policy.R"), local = TRUE)

make_opt <- function(locked = FALSE,
                     enable_ai = TRUE,
                     providers = c("bigomics", "openai", "anthropic", "google",
                                   "github", "mistral", "custom"),
                     policy = .opg_ai_read_policy(file.path(.repo_dir, "etc/ai_model_policy.json"))) {
  models <- .opg_ai_build_models(policy, providers)
  list(
    ENABLE_AI                = enable_ai,
    AI_PROVIDER_LOCKED       = locked,
    AI_PROVIDERS             = providers,
    AI_MENU_REPORTS          = .opg_ai_menu_allowlist(models, providers, "reports"),
    AI_MENU_IMAGES           = .opg_ai_menu_allowlist(models, providers, "images"),
    AI_MENU_COPILOT_DEEP     = .opg_ai_menu_allowlist(models, providers, "copilot_deep"),
    AI_MENU_COPILOT_BALANCED = .opg_ai_menu_allowlist(models, providers, "copilot_balanced"),
    AI_MODELS                = models
  )
}
opt <- make_opt()

.board_dir <- if (dir.exists("components/board.user/R")) {
  "components/board.user/R"
} else {
  "../../components/board.user/R"
}
source(file.path(.board_dir, "appsettings_server.R"), local = TRUE)

make_auth <- function(admin = FALSE) {
  shiny::reactiveValues(logged = TRUE, ADMIN = admin, user_dir = tempdir())
}

# ===========================================================================
# credential store
# ===========================================================================

test_that("a non-bigomics provider + key stores a nullary closure returning the key", {
  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(), pgx = shiny::reactiveValues()), {
      session$setInputs(ai_provider = "openai", ai_api_key = "sk-secret-123")
      cred <- session$userData[["ai_credentials"]]
      expect_true(is.function(cred))
      expect_equal(length(formals(cred)), 0L)
      expect_equal(cred(), "sk-secret-123")
      expect_equal(session$userData[["ai_provider"]], "openai")
    }
  )
})

test_that("a non-bigomics provider with an empty key stores no closure", {
  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(), pgx = shiny::reactiveValues()), {
      session$setInputs(ai_provider = "openai", ai_api_key = "")
      expect_null(session$userData[["ai_credentials"]])
    }
  )
})

test_that("provider = bigomics stores NULL credentials (env-var path preserved)", {
  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(), pgx = shiny::reactiveValues()), {
      session$setInputs(ai_provider = "openai", ai_api_key = "sk-x")
      expect_true(is.function(session$userData[["ai_credentials"]]))
      session$setInputs(ai_provider = "bigomics")
      expect_null(session$userData[["ai_credentials"]])
      expect_equal(session$userData[["ai_provider"]], "bigomics")
    }
  )
})

test_that("the custom provider stores its base_url alongside the credential", {
  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(), pgx = shiny::reactiveValues()), {
      session$setInputs(
        ai_provider = "custom", ai_api_key = "ck", ai_base_url = "https://ep/v1"
      )
      expect_equal(session$userData[["ai_credentials"]](), "ck")
      expect_equal(session$userData[["ai_base_url"]], "https://ep/v1")
    }
  )
})

# ===========================================================================
# legacy option compatibility
# ===========================================================================

test_that("model selections still write the legacy llm_model / img_model options", {
  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(), pgx = shiny::reactiveValues()), {
      session$setInputs(
        enable_ai = TRUE,
        llm_reports = "gpt-5.4-nano", llm_images = "dall-e-3",
        llm_copilot_deep = "gpt-5.4-mini",
        llm_copilot_balanced = "gpt-5.4-nano"
      )
      expect_equal(session$userData[["llm_model"]], "gpt-5.4-nano")
      expect_equal(session$userData[["img_model"]], "dall-e-3")
      expect_equal(session$userData[["llm_copilot_deep"]], "gpt-5.4-mini")
      expect_equal(session$userData[["llm_copilot_balanced"]], "gpt-5.4-nano")
    }
  )
})

# ===========================================================================
# provider-aware menu repopulation
# ===========================================================================

test_that("OPG policy builds provider menus from the omicsai catalog APIs", {
  expect_true("github" %in% opt$AI_PROVIDERS)
  expect_true("mistral" %in% opt$AI_PROVIDERS)
  expect_equal(opt$AI_MODELS$github$reports,
               c("openai/gpt-4.1", "openai/gpt-4o", "openai/gpt-4o-mini"))
  expect_true("mistral-medium-latest" %in% opt$AI_MODELS$mistral$reports)
  expect_equal(opt$AI_MODELS$mistral$copilot_balanced[[1]],
               "mistral-medium-latest")
})

test_that("OPG policy can disable a provider and reorder menu defaults", {
  policy <- .opg_ai_read_policy(file.path(.repo_dir, "etc/ai_model_policy.json"))
  policy$providers$mistral$menus$reports$prefer <- "mistral-small-latest"
  limited <- make_opt(providers = c("openai", "mistral"), policy = policy)

  expect_equal(limited$AI_PROVIDERS, c("openai", "mistral"))
  expect_equal(limited$AI_MODELS$mistral$reports[[1]],
               "mistral-small-latest")
  expect_null(make_opt(providers = "openai")$AI_MODELS$mistral)
})

test_that("BigOmics defaults stay unchanged and model menus remain hidden", {
  expect_equal(opt$AI_MODELS$bigomics$reports[[1]], "openai:gpt-5.4-nano")
  expect_equal(opt$AI_MODELS$bigomics$images[[1]],
               "gemini-3.1-flash-image-preview")
  expect_equal(opt$AI_MODELS$bigomics$copilot_deep,
               "openai:gpt-5.4-mini")
  expect_equal(opt$AI_MODELS$bigomics$copilot_balanced,
               "openai:gpt-5.4-nano")

  ui_source <- readLines(file.path(.repo_dir,
                                   "components/board.user/R/appsettings_ui.R"))
  expect_true(any(grepl("input.ai_provider != 'bigomics'", ui_source,
                        fixed = TRUE)))
  expect_true(any(grepl("ai_test_status", ui_source, fixed = TRUE)))
})

test_that("changing provider repopulates the four menus from that provider's catalog", {
  captured <- list()
  testthat::local_mocked_bindings(
    updateSelectInput = function(session, inputId, ...) {
      captured[[inputId]] <<- list(
        choices = list(...)$choices,
        selected = list(...)$selected
      )
    },
    .package = "shiny"
  )
  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(), pgx = shiny::reactiveValues()), {
      session$setInputs(ai_provider = "openai")
      expect_equal(captured$llm_reports$choices,
                   c("gpt-5.4-nano", "gpt-5.4-mini",
                     "gpt-4o", "gpt-4o-mini"))
      expect_equal(captured$llm_reports$selected, "gpt-5.4-nano")
      expect_equal(captured$llm_images$choices, c("dall-e-3", "dall-e-2"))
      expect_equal(captured$llm_copilot_deep$selected, "gpt-5.4-mini")
      expect_equal(captured$llm_copilot_balanced$selected, "gpt-5.4-nano")
    }
  )
})

test_that("Test & load narrows menus to live models allowed by OPG policy", {
  captured <- list()
  requested <- list()
  alert <- NULL
  testthat::local_mocked_bindings(
    updateSelectInput = function(session, inputId, ...) {
      captured[[inputId]] <<- list(
        choices = list(...)$choices,
        selected = list(...)$selected
      )
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    shinyalert = function(...) {
      alert <<- list(...)
    },
    .package = "shinyalert"
  )
  testthat::local_mocked_bindings(
    ai.list_provider_models = function(provider, key, base_url = NULL) {
      requested <<- list(provider = provider, key = key, base_url = base_url)
      c("gpt-4o-mini", "not-in-opg-policy")
    },
    .package = "omicsai"
  )

  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(), pgx = shiny::reactiveValues()), {
      session$setInputs(ai_provider = "openai", ai_api_key = "sk-live")
      session$setInputs(ai_test_load = 1)

      expect_equal(requested$provider, "openai")
      expect_equal(requested$key, "sk-live")
      expect_equal(captured$llm_reports$choices, "gpt-4o-mini")
      expect_equal(captured$llm_reports$selected, "gpt-4o-mini")
      expect_equal(captured$llm_copilot_deep$choices, "gpt-4o-mini")
      expect_equal(captured$llm_images$choices, c("dall-e-3", "dall-e-2"))
      expect_equal(alert$title, "API key is correctly set")
      expect_equal(alert$type, "success")
      expect_equal(ai_test_status()$state, "ok")
      expect_equal(ai_test_status()$label, "OK")
    }
  )
})

test_that("Test & load falls back to static catalog choices on empty live results", {
  captured <- list()
  alert <- NULL
  testthat::local_mocked_bindings(
    updateSelectInput = function(session, inputId, ...) {
      captured[[inputId]] <<- list(...)$choices
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    shinyalert = function(...) {
      alert <<- list(...)
    },
    .package = "shinyalert"
  )
  testthat::local_mocked_bindings(
    ai.list_provider_models = function(...) character(0),
    .package = "omicsai"
  )

  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(), pgx = shiny::reactiveValues()), {
      session$setInputs(ai_provider = "openai", ai_api_key = "sk-live")
      session$setInputs(ai_test_load = 1)

      expect_equal(captured$llm_reports,
                   c("gpt-5.4-nano", "gpt-5.4-mini",
                     "gpt-4o", "gpt-4o-mini"))
      expect_equal(captured$llm_images, c("dall-e-3", "dall-e-2"))
      expect_equal(alert$title, "Could not load provider models")
      expect_equal(alert$type, "error")
      expect_equal(ai_test_status()$state, "error")
      expect_equal(ai_test_status()$label, "Error")
    }
  )
})

test_that("Test & load passes the custom provider base_url and loads live text menus", {
  requested <- NULL
  captured <- list()
  alert <- NULL
  testthat::local_mocked_bindings(
    updateSelectInput = function(session, inputId, ...) {
      captured[[inputId]] <<- list(...)$choices
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    shinyalert = function(...) {
      alert <<- list(...)
    },
    .package = "shinyalert"
  )
  testthat::local_mocked_bindings(
    ai.list_provider_models = function(provider, key, base_url = NULL) {
      requested <<- list(provider = provider, key = key, base_url = base_url)
      c("custom-live-model", "custom-large")
    },
    .package = "omicsai"
  )

  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(), pgx = shiny::reactiveValues()), {
      session$setInputs(
        ai_provider = "custom",
        ai_api_key = "ck",
        ai_base_url = "https://llm.example/v1"
      )
      session$setInputs(ai_test_load = 1)

      expect_equal(requested$provider, "custom")
      expect_equal(requested$key, "ck")
      expect_equal(requested$base_url, "https://llm.example/v1")
      expect_equal(captured$llm_reports,
                   c("custom-live-model", "custom-large"))
      expect_equal(captured$llm_copilot_deep,
                   c("custom-live-model", "custom-large"))
      expect_equal(captured$llm_copilot_balanced,
                   c("custom-live-model", "custom-large"))
      expect_equal(captured$llm_images, character(0))
      expect_equal(alert$type, "success")
      expect_equal(ai_test_status()$state, "ok")
    }
  )
})

test_that("a provider whose menu is empty repopulates with no choices", {
  captured <- list()
  testthat::local_mocked_bindings(
    updateSelectInput = function(session, inputId, ...) {
      captured[[inputId]] <<- list(...)$choices
    },
    .package = "shiny"
  )
  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(), pgx = shiny::reactiveValues()), {
      session$setInputs(ai_provider = "anthropic")
      expect_equal(captured$llm_images, character(0))
      expect_equal(captured$llm_reports,
                   c("claude-opus-4-8", "claude-sonnet-4-6",
                     "claude-haiku-4-5"))
    }
  )
})

test_that("GitHub and Mistral providers repopulate from their catalogs", {
  captured <- list()
  testthat::local_mocked_bindings(
    updateSelectInput = function(session, inputId, ...) {
      captured[[inputId]] <<- list(...)$choices
    },
    .package = "shiny"
  )
  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(), pgx = shiny::reactiveValues()), {
      session$setInputs(ai_provider = "github")
      expect_equal(captured$llm_reports,
                   c("openai/gpt-4.1", "openai/gpt-4o", "openai/gpt-4o-mini"))
      expect_equal(captured$llm_images, character(0))
      expect_equal(captured$llm_copilot_deep,
                   c("openai/gpt-4.1", "openai/gpt-4o", "openai/gpt-4o-mini"))
      expect_equal(captured$llm_copilot_balanced,
                   c("openai/gpt-4o-mini", "openai/gpt-4.1", "openai/gpt-4o"))

      session$setInputs(ai_provider = "mistral")
      expect_equal(captured$llm_reports,
                   c("mistral-large-latest", "mistral-medium-latest",
                     "mistral-small-latest", "ministral-8b-latest"))
      expect_equal(captured$llm_images, character(0))
      expect_equal(captured$llm_copilot_deep,
                   c("mistral-large-latest", "mistral-medium-latest",
                     "mistral-small-latest", "ministral-8b-latest"))
      expect_equal(captured$llm_copilot_balanced,
                   c("mistral-medium-latest", "mistral-large-latest",
                     "mistral-small-latest", "ministral-8b-latest"))
    }
  )
})

test_that("a provider missing from the menu policy falls back to the union allowlist", {
  captured <- list()
  testthat::local_mocked_bindings(
    updateSelectInput = function(session, inputId, ...) {
      captured[[inputId]] <<- list(...)$choices
    },
    .package = "shiny"
  )
  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(), pgx = shiny::reactiveValues()), {
      session$setInputs(ai_provider = "unknown")
      expect_equal(captured$llm_reports, opt$AI_MENU_REPORTS)
      expect_equal(captured$llm_images, opt$AI_MENU_IMAGES)
    }
  )
})

# ===========================================================================
# admin lock
# ===========================================================================

test_that("a non-admin sees a disabled provider dropdown when AI_PROVIDER_LOCKED", {
  opt <<- make_opt(locked = TRUE)
  on.exit(opt <<- make_opt(), add = TRUE)
  toggled <- list(disabled = character(0), enabled = character(0))
  testthat::local_mocked_bindings(
    disable = function(id, ...) toggled$disabled <<- c(toggled$disabled, id),
    enable  = function(id, ...) toggled$enabled  <<- c(toggled$enabled, id),
    .package = "shinyjs"
  )
  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(admin = FALSE), pgx = shiny::reactiveValues()), {
      session$flushReact()
      expect_true("ai_provider" %in% toggled$disabled)
      expect_false("ai_provider" %in% toggled$enabled)
    }
  )
})

test_that("an admin keeps the provider dropdown enabled even when locked", {
  opt <<- make_opt(locked = TRUE)
  on.exit(opt <<- make_opt(), add = TRUE)
  toggled <- list(disabled = character(0), enabled = character(0))
  testthat::local_mocked_bindings(
    disable = function(id, ...) toggled$disabled <<- c(toggled$disabled, id),
    enable  = function(id, ...) toggled$enabled  <<- c(toggled$enabled, id),
    .package = "shinyjs"
  )
  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(admin = TRUE), pgx = shiny::reactiveValues()), {
      session$flushReact()
      expect_true("ai_provider" %in% toggled$enabled)
      expect_false("ai_provider" %in% toggled$disabled)
    }
  )
})

# ===========================================================================
# enable-AI gate: publishes ai_enabled + toggles the AI tabs
# ===========================================================================

test_that("enabling AI publishes ai_enabled=TRUE and exposes enable_ai=TRUE", {
  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(), pgx = shiny::reactiveValues()), {
      session$setInputs(enable_ai = TRUE)
      session$flushReact()
      expect_true(isTRUE(session$userData[["ai_enabled"]]))
      expect_true(isTRUE(session$returned$enable_ai()))
    }
  )
})

test_that("disabling the AI switch publishes ai_enabled=FALSE and exposes enable_ai=FALSE", {
  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(), pgx = shiny::reactiveValues()), {
      session$setInputs(enable_ai = FALSE)
      session$flushReact()
      expect_false(isTRUE(session$userData[["ai_enabled"]]))
      expect_false(isTRUE(session$returned$enable_ai()))
    }
  )
})

test_that("an unlicensed deployment (ENABLE_AI=FALSE) forces ai_enabled off even when the switch is on", {
  opt <<- make_opt(enable_ai = FALSE)
  on.exit(opt <<- make_opt(), add = TRUE)
  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(), pgx = shiny::reactiveValues()), {
      session$setInputs(enable_ai = TRUE)
      session$flushReact()
      expect_false(isTRUE(session$userData[["ai_enabled"]]))
    }
  )
})

# ===========================================================================
# logout clears the session credential
# ===========================================================================

test_that("logging out clears the stored credential and resets the provider", {
  shiny::testServer(AppSettingsBoard,
    args = list(id = "s", auth = make_auth(), pgx = shiny::reactiveValues()), {
      session$setInputs(ai_provider = "openai", ai_api_key = "sk-leak")
      expect_true(is.function(session$userData[["ai_credentials"]]))

      auth$logged <- FALSE
      session$flushReact()

      expect_null(session$userData[["ai_credentials"]])
      expect_equal(session$userData[["ai_provider"]], "bigomics")
    }
  )
})
