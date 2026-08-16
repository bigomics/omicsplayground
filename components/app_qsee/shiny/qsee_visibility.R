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
  purge_on_hide <- qsee_opt_enabled("ENABLE_PLOTLY_PURGE", TRUE)
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
        var purgeOnHide = %s;
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
        function purgeNearby() {
          if (!purgeOnHide) return;
          var root = el.closest('.big-tab') || el.parentElement;
          if (!root) return;
          var nodes = root.querySelectorAll(
            '.js-plotly-plot, .plotly.html-widget, .iheatmapr'
          );
          for (var i = 0; i < nodes.length; i++) {
            var node = nodes[i];
            if (!node.querySelector('svg, canvas')) continue;
            if (window.Plotly && typeof Plotly.purge === 'function') {
              try { Plotly.purge(node); } catch (e) {}
            }
            while (node.firstChild) node.removeChild(node.firstChild);
          }
        }
        function check() {
          // offsetParent is null while a parent nav panel is display:none.
          var visible = el.offsetParent !== null && el.getClientRects().length > 0;
          if (!visible) purgeNearby();
          report(visible);
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
      ns("visible_probe"),
      ns("is_visible"),
      if (isTRUE(purge_on_hide)) "true" else "false"
    )))
  )
}

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

#' Purge this board's Plotly/iheatmapr widgets while hidden; redraw on return.
#'
#' Drops the drawn SVG/WebGL tree so hidden Qsee tabs do not keep tens of
#' thousands of DOM nodes. Input widgets stay mounted, so Shiny input state
#' is preserved. Returning to the tab increments a tick that plot renderers
#' must read (via [qsee_with_redraw()] or a direct `tick()` call).
#'
#' Toggle with `ENABLE_PLOTLY_PURGE` in etc/OPTIONS. When disabled, the
#' returned reactive is a constant and never invalidates.
#'
#' @param is_visible reactive logical from [qsee_is_visible()]
#' @param session module `session` (used to namespace the purge selector)
#' @param label optional debug label
#' @return a `reactive` integer
qsee_plotly_purge <- function(is_visible, session, label = NULL) {
  redraw_tick <- shiny::reactiveVal(0L)

  if (!qsee_opt_enabled("ENABLE_PLOTLY_PURGE", TRUE)) {
    if (!is.null(label)) {
      message("[", label, "] plotly purge DISABLED")
    }
    return(shiny::reactive(redraw_tick()))
  }

  ## Only redraw after a real hide. The first visible=TRUE (first visit
  ## or the initially selected tab) should not force a second render.
  ever_hidden <- FALSE

  shiny::observeEvent(
    is_visible(),
    {
      visible <- isTRUE(is_visible())
      if (!visible) {
        ever_hidden <<- TRUE
        prefix <- session$ns("")
        shinyjs::runjs(sprintf(
          "(function() {
            var prefix = '%s';
            var nodes = document.querySelectorAll(
              '.js-plotly-plot, .plotly.html-widget, .iheatmapr'
            );
            for (var i = 0; i < nodes.length; i++) {
              var el = nodes[i];
              if (prefix && el.id && el.id.indexOf(prefix) !== 0) continue;
              if (window.Plotly && typeof Plotly.purge === 'function') {
                try { Plotly.purge(el); } catch (e) {}
              }
              while (el.firstChild) el.removeChild(el.firstChild);
            }
          })();",
          prefix
        ))
        if (!is.null(label)) {
          message("[", label, "] plotly purged")
        }
      } else if (ever_hidden) {
        redraw_tick(shiny::isolate(redraw_tick()) + 1L)
        if (!is.null(label)) {
          message("[", label, "] plotly redraw requested")
        }
      }
    },
    ignoreInit = TRUE
  )

  shiny::reactive(redraw_tick())
}

#' Keep the bigLoaders spinner up until Plotly has actually painted.
#'
#' `PlotModuleUI()` wraps every plot in `bigLoaders::useSpinner()`, whose
#' spinner.js hides the spinner on `shiny:value`. That event fires when the
#' figure data reaches the browser, *before* plotly.js draws it -- so the card
#' sits blank for the whole client-side draw, which for a big scatter or
#' heatmap is the slow part. Same story after [qsee_plotly_purge()] empties a
#' plot: the redraw is client-side only.
#'
#' This re-shows the spinner at `shiny:value` and hides it once the output
#' element actually contains a drawn `svg`/`canvas`. Include once, in the
#' module UI. Scoped by namespace prefix so it cannot affect other boards.
#'
#' Only plotly/iheatmapr outputs are touched. `plotlib = "base"` plots are
#' rendered to an image server-side and have no client-side draw phase.
#'
#' @param ns module namespace function (`shiny::NS(id)`)
#' @return a `shiny.tag` script element
qsee_plot_spinner_js <- function(ns) {
  shiny::tags$script(shiny::HTML(sprintf(
    "(function() {
      var prefix = '%s';
      window.__qseeSpinnerPrefixes = window.__qseeSpinnerPrefixes || {};
      if (window.__qseeSpinnerPrefixes[prefix]) return;
      window.__qseeSpinnerPrefixes[prefix] = true;

      // ~10s at 60fps. A plot that never paints (empty data, render error)
      // must not leave the spinner spinning forever.
      var MAX_FRAMES = 600;

      function isWidget(el) {
        return el.classList.contains('plotly') ||
               el.classList.contains('js-plotly-plot') ||
               el.classList.contains('iheatmapr');
      }
      function painted(el) {
        return !!el.querySelector('svg, canvas');
      }

      $(document).on('shiny:value', function(e) {
        if (!e.name || e.name.indexOf(prefix) !== 0) return;
        var el = document.getElementById(e.name);
        if (!el || !isWidget(el)) return;
        var sp = document.getElementById(e.name + '-spinner');
        if (!sp) return;
        if (painted(el)) return;

        // Runs after spinner.js's own handler (its script is a page-level
        // htmlDependency, registered first), so this undoes its early hide.
        $(sp).show();
        // visibility, not display: plotly measures the container to size the
        // figure, so it has to keep its layout box while drawing.
        el.style.visibility = 'hidden';

        var frames = 0;
        (function wait() {
          if (painted(el) || ++frames > MAX_FRAMES) {
            $(sp).hide();
            el.style.visibility = 'inherit';
            return;
          }
          window.requestAnimationFrame(wait);
        })();
      });
    })();",
    ns("")
  )))
}

#' Make a plot renderer depend on a [qsee_plotly_purge()] tick.
qsee_with_redraw <- function(tick, func) {
  force(func)
  force(tick)
  function(...) {
    tick()
    func(...)
  }
}

## NOTE: boards deliberately have no lazy-mount helper. Each board UI puts
## its body behind `shiny::uiOutput(ns("lazy_body"))` and the server fills it
## with a plain `renderUI()`. Shiny suspends hidden outputs by default
## (`suspendWhenHidden = TRUE`), and a board sitting in an inactive bigdash
## `.big-tab` is `display:none`, so that renderUI does not run until the tab
## is first shown -- and the built body then stays in the DOM. Pair with
## [qsee_plotly_purge()] to drop the drawn Plotly trees while hidden.

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
