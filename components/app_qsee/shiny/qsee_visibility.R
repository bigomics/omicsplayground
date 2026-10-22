##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Tab-visibility gating and Plotly purging for QSEE boards is provided by
## bigdash's visibility toolbox (bigdash::bd_visibility_probe(),
## bigdash::bd_is_visible(), bigdash::bd_redraw_tick(),
## bigdash::bd_with_redraw()) -- see R/visibility.R in the bigdash package.
##
## What is left here is the one bit that is app-specific: whether purging is
## enabled at all, driven by ENABLE_PLOTLY_PURGE in etc/OPTIONS.
##

#' Read a boolean from the global `opt` list (etc/OPTIONS) or `getOption()`.
qsee_opt_enabled <- function(name, default = TRUE) {
  if (exists("opt", envir = .GlobalEnv, inherits = FALSE)) {
    opt_list <- get("opt", envir = .GlobalEnv, inherits = FALSE)
    if (is.list(opt_list) && !is.null(opt_list[[name]])) {
      return(isTRUE(opt_list[[name]]))
    }
  }
  isTRUE(getOption(name, default))
}

#' Whether Plotly purge-while-hidden is enabled (`ENABLE_PLOTLY_PURGE`).
#'
#' Pass straight through as `bigdash::bd_is_visible(input, purge = ...)`'s
#' `purge` argument.
qsee_purge_enabled <- function() {
  qsee_opt_enabled("ENABLE_PLOTLY_PURGE", TRUE)
}

#' Resolve a board's `purge` setting: an explicit override, or the option.
#'
#' `qsee_server(purge = ...)` threads a single value down through every
#' board so it can be flipped for the whole app at launch (e.g. for an
#' on/off A-B comparison from the test harness) without touching
#' etc/OPTIONS. `NULL` (the default all the way down) keeps the existing
#' per-install [qsee_purge_enabled()] behaviour; `TRUE`/`FALSE` overrides it.
#'
#' @param purge `NULL`, or an explicit `TRUE`/`FALSE` override.
qsee_resolve_purge <- function(purge = NULL) {
  if (is.null(purge)) qsee_purge_enabled() else isTRUE(purge)
}

## NOTE: boards deliberately have no lazy-mount helper. Each board UI puts
## its body behind `shiny::uiOutput(ns("lazy_body"))` and the server fills it
## with a plain `renderUI()`. Shiny suspends hidden outputs by default
## (`suspendWhenHidden = TRUE`), and a board sitting in an inactive bigdash
## `.big-tab` is `display:none`, so that renderUI does not run until the tab
## is first shown -- and the built body then stays in the DOM. Pair with
## `bigdash::bd_is_visible(input, purge = qsee_resolve_purge(purge))` to drop
## the drawn Plotly trees while hidden.

## NOTE: boards also have no precomputation-cache helper. A board's heavy
## work is a plain `shiny::reactive()`, which already caches (returning to a
## tab reuses the last result), invalidates when its inputs change, and does
## nothing while the board is hidden -- its only consumers are plot/table
## outputs, which Shiny suspends. An input change that happens while hidden
## invalidates it but is only computed on the next visit.
##
## Where the body must read a reactive *without* depending on it (e.g. the
## batch-effects board reads input$main_param but only recomputes on the
## Recompute button), use `shiny::eventReactive()` -- that is exactly a
## dependency expression plus an isolated body.
##
## This all holds only as long as every consumer is lazy. Reading such a
## reactive from an `observe()`/`observeEvent()`, or from an output with
## `suspendWhenHidden = FALSE`, pulls it eagerly and the board computes while
## hidden. Where that is unavoidable, guard the observer with
## `shiny::req(is_visible())` -- see qsee_batchcorrect_server.R.
