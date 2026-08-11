##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Shared tab-visibility helpers for QSEE boards. Each board UI embeds
## qsee_visibility_probe(); the matching server reads input$is_visible
## (via qsee_is_visible) and gates expensive precomputation on it.
##

#' Invisible IntersectionObserver probe for board tab visibility.
#'
#' Reports whether this board's UI is actually shown on screen (e.g. its
#' sidebar tab is selected) via `input$is_visible`, so the server can
#' pause/resume expensive reactives. Works generically via
#' IntersectionObserver instead of depending on any particular
#' tab-container implementation (bigdash, bslib navset, etc.).
#'
#' NOTE: the probe must stay in normal layout flow (no display:none) or
#' it can never register as visible. It is transparent rather than hidden so
#' both IntersectionObserver and the fallback layout check can see it.
#'
#' @param ns module namespace function (`shiny::NS(id)`)
#' @return a `shiny.tag.list` to prepend to the board UI
qsee_visibility_probe <- function(ns) {
  shiny::tagList(
    shiny::tags$div(
      id = ns("visible_probe"),
      style = "width:1px; height:1px; opacity:0; pointer-events:none;"
    ),
    shiny::tags$script(shiny::HTML(sprintf(
      "(function() {
        var el = document.getElementById('%s');
        if (!el) return;
        var inputId = '%s';
        var lastState = null;
        function report(state) {
          // Shiny.setInputValue is not attached until shiny.js finishes
          // initializing, which can be after this script runs (it executes
          // as soon as its containing markup is parsed/inserted). Retrying
          // is safe: on failure we leave lastState untouched so the very
          // same state is reported again on the next check() -- including
          // the guaranteed one on 'shiny:connected' below -- instead of
          // being silently dropped forever.
          if (!(window.Shiny && typeof Shiny.setInputValue === 'function'))
            return;
          if (state !== lastState) {
            lastState = state;
            Shiny.setInputValue(inputId, state, {priority: 'event'});
          }
        }
        function check() {
          // offsetParent is null while a parent nav panel is display:none.
          report(el.offsetParent !== null && el.getClientRects().length > 0);
        }
        if (window.IntersectionObserver) {
          var io = new IntersectionObserver(function() { check(); });
          io.observe(el);
        }
        // bslib changes classes/attributes on its parent nav panel. An
        // initially hidden child does not reliably receive an intersection
        // callback when that panel is activated, so observe those changes too.
        var mo = new MutationObserver(function() {
          window.requestAnimationFrame(check);
        });
        mo.observe(document.documentElement, {
          attributes: true,
          attributeFilter: ['class', 'style', 'hidden', 'aria-hidden'],
          childList: true,
          subtree: true
        });
        check();
        document.addEventListener('shiny:connected', check);
      })();",
      ns("visible_probe"), ns("is_visible")
    )))
  )
}

#' Reactive that is TRUE only while this board's tab is on screen.
#'
#' @param input module `input`
#' @param label optional label for debug messages (e.g. module name)
#' @return a `reactive` logical
qsee_is_visible <- function(input, label = NULL) {
  shiny::reactive({
    if (!is.null(label)) {
      message("[", label, "] visible = ", input$is_visible)
    }
    isTRUE(input$is_visible)
  })
}

#' Visibility-gated board cache: recompute only when inputs change.
#'
#' Visibility controls *when* work is allowed to run, but does **not**
#' invalidate the cache. Returning to a nav tab reuses the last result
#' unless `deps` changed (including changes that happened while hidden —
#' those are applied on the next visit).
#'
#' @param is_visible reactive logical from [qsee_is_visible()]
#' @param deps zero-arg function or reactive that returns a value whose
#'   change marks the cache stale (e.g. `function() list(rX(), rY())`)
#' @param compute zero-arg function that performs the heavy work and
#'   returns the result object. Evaluated under `isolate()` so it does
#'   not create extra reactive dependencies.
#' @param label optional debug label
#' @return a `reactive` yielding the latest computed result
qsee_board_cache <- function(is_visible, deps, compute, label = NULL) {
  force(compute)
  if (!is.function(deps)) {
    stop("qsee_board_cache: `deps` must be a zero-arg function or reactive")
  }

  result <- shiny::reactiveVal(NULL)
  stale <- shiny::reactiveVal(TRUE)
  ## A data set can arrive while the cache is already stale (for example,
  ## when the app starts without a PGX).  Toggling `stale` to TRUE again
  ## would not invalidate a reactive observer, so retain a monotonically
  ## changing generation as an explicit invalidation signal.
  generation <- shiny::reactiveVal(0L)

  shiny::observeEvent(
    deps(),
    {
      stale(TRUE)
      generation(shiny::isolate(generation()) + 1L)
      if (!is.null(label)) {
        message("[", label, "] inputs changed — cache marked stale")
      }
    },
    ignoreNULL = TRUE
  )

  shiny::observe({
    generation()
    shiny::req(is_visible())
    shiny::req(isTRUE(stale()))
    if (!is.null(label)) {
      message("[", label, "] computing...")
    }
    ## isolate so input reads inside compute don't re-trigger this observe
    res <- shiny::isolate(compute())
    result(res)
    stale(FALSE)
  })

  shiny::reactive({
    shiny::req(!is.null(result()))
    result()
  })
}
