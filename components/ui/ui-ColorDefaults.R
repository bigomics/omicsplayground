##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## ---------------------------------------------------------------
## Color Theme Defaults
## ---------------------------------------------------------------

COLOR_THEME_DEFAULTS <- list(
  primary    = "#f23451", # Up / High
  secondary  = "#3181de", # Down / Low / scatter_color
  neutral    = "#eeeeee", # Mid (heatmap midpoint)
  ns_color   = "#eeeeee", # Not significant dots
  bar_color  = "#A6CEE3", # Bar color
  accent     = "#e3a45a", # Significant in one (color_one)
  success    = "#5B9B5B", # Significant in both (color_both)
  line       = "#00EE00", # Enrichment line
  palette    = "default",
  palette_c1 = "#3181de", # Custom gradient start
  palette_c2 = "#eeeeee", # Custom gradient middle
  palette_c3 = "#f23451" # Custom gradient end
)

## Maps theme keys to the editor input IDs they control
COLOR_THEME_MAPPING <- list(
  primary   = c("color_up", "color_high"),
  secondary = c("color_down", "color_low", "scatter_color", "rank_color_line"),
  neutral   = c("color_mid"),
  ns_color  = c("color_ns"),
  bar_color = c("bar_color"),
  accent    = c("color_one"),
  success   = c("color_both"),
  line      = c("color_line")
)

## ---------------------------------------------------------------
## Per-session reactive store
## ---------------------------------------------------------------
##
## The store must NOT be a process-level singleton. reactiveValues() registers
## an onDestroy handler on whatever reactive domain is current when it is
## created, and any later write raises "Can't access reactive ...; its module
## session has been destroyed". A store built during the first browser session
## and cached for the process is therefore already dead for the second one:
## reload the page and the first write (appsettings_server.R, on auth$logged)
## fails. It also leaked one user's theme into the next session.
##
## So: one store per browser session, kept on the session's userData.

.color_theme_env <- new.env(parent = emptyenv())
.color_theme_env$theme <- NULL ## fallback for use outside any session

.new_color_theme <- function() {
  do.call(shiny::reactiveValues, COLOR_THEME_DEFAULTS)
}

#' The colour theme store for the current browser session.
#'
#' Anchored on the *root* session rather than the calling module session: the
#' theme is app-wide, and a module session can be torn down (its UI removed)
#' long before the browser session ends, which would kill the store early.
#' Module sessions share the root session's `userData`, so every board sees
#' the same object.
#'
#' @param session reactive domain to anchor on; defaults to the current one
#' @return reactiveValues
color_theme_store <- function(session = shiny::getDefaultReactiveDomain()) {
  if (is.null(session)) {
    ## No reactive domain: console, tests, or UI built outside a session. A
    ## reactiveValues created with no domain registers no onDestroy handler,
    ## so it is never marked destroyed and is safe to keep for the process.
    if (is.null(.color_theme_env$theme)) {
      .color_theme_env$theme <- .new_color_theme()
    }
    return(.color_theme_env$theme)
  }

  root <- if (is.function(session$rootScope)) session$rootScope() else session
  if (is.null(root$userData$color_theme)) {
    root$userData$color_theme <- shiny::withReactiveDomain(root, .new_color_theme())
  }
  root$userData$color_theme
}

#' Initialise the colour theme reactive store (call once in server).
#' @return invisible reactiveValues
init_color_theme <- function() {
  invisible(color_theme_store())
}

#' Return the colour theme reactive store.
#' @return reactiveValues
get_color_theme <- function() {
  color_theme_store()
}

#' Save the colour theme to a user directory as JSON.
save_color_theme <- function(theme_list, user_dir) {
  if (is.null(user_dir) || !dir.exists(user_dir)) {
    return(invisible(NULL))
  }
  path <- file.path(user_dir, "color_theme.json")
  jsonlite::write_json(theme_list, path, auto_unbox = TRUE, pretty = TRUE)
  invisible(path)
}

#' Load colour theme from a user directory. Returns NULL if file missing/invalid.
load_color_theme <- function(user_dir) {
  if (is.null(user_dir)) {
    return(NULL)
  }
  path <- file.path(user_dir, "color_theme.json")
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(
    jsonlite::read_json(path, simplifyVector = TRUE),
    error = function(e) NULL
  )
}

#' Push theme colour changes into the plot editor's colour inputs.
#'
#' Registered as the `bigdash.editor_theme_observer` hook (see
#' \code{bigdash::bd_hook}) and called by PlotModuleServer with the calling
#' module's session.  Lives here rather than in bigdash because
#' COLOR_THEME_MAPPING names the editor's own input IDs.
#'
#' @param parent_session Session of the board module owning the editor inputs.
plotmodule_theme_observer <- function(parent_session) {
  reg <- .theme_target_registry(parent_session)

  ## Register this module, then push the current theme to it alone. This
  ## replaces the old `ignoreInit = FALSE`, which existed so a lazily loaded
  ## module would pick up values changed before it was navigated to.
  reg$targets[[length(reg$targets) + 1L]] <- parent_session
  .push_theme_to(list(parent_session), shiny::isolate(
    shiny::reactiveValuesToList(get_color_theme())
  ))

  ## The observers themselves are installed once per browser session, not once
  ## per plot module. Previously every module installed its own set: with 70
  ## modules passing parent_session and 9 observers each (8 colour keys plus
  ## the palette), that was ~630 observers, all firing on materialisation and
  ## issuing a separate runjs() per mapped input id.
  if (isTRUE(reg$installed)) {
    return(invisible(NULL))
  }
  reg$installed <- TRUE

  theme <- get_color_theme()

  lapply(names(COLOR_THEME_MAPPING), function(key) {
    shiny::observeEvent(theme[[key]],
      {
        vals <- stats::setNames(list(theme[[key]]), key)
        .push_theme_to(reg$targets, vals)
      },
      ignoreInit = TRUE ## new modules are seeded at registration, above
    )
  })

  shiny::observeEvent(theme$palette,
    {
      for (s in reg$targets) {
        shiny::updateSelectInput(s, "palette", selected = theme$palette)
      }
    },
    ignoreInit = TRUE
  )

  invisible(NULL)
}

#' Per-session registry of plot modules wanting theme updates.
#'
#' An environment, so the observers installed on the first registration keep
#' seeing modules that register later.
.theme_target_registry <- function(session) {
  root <- if (is.function(session$rootScope)) session$rootScope() else session
  if (is.null(root$userData$theme_targets)) {
    reg <- new.env(parent = emptyenv())
    reg$targets <- list()
    reg$installed <- FALSE
    root$userData$theme_targets <- reg
  }
  root$userData$theme_targets
}

#' Push theme values to plot modules.
#'
#' `Shiny.setInputValue` is used rather than only `updateColourInput` because
#' the colourpicker widget may not be initialised while the editor modal is
#' hidden; the updateColourInput call keeps the swatch right for when it is
#' opened.
#'
#' Every setInputValue for this push goes out in ONE runjs() call. The old code
#' issued one per module per input id -- on the order of 840 calls as boards
#' materialised, each landing as its own message and its own input change, so
#' each one could trigger another reactive flush.
#'
#' @param targets list of module sessions
#' @param vals named list of theme key -> colour
.push_theme_to <- function(targets, vals) {
  if (!length(targets) || !length(vals)) {
    return(invisible(NULL))
  }
  js <- character(0)
  for (key in names(vals)) {
    input_ids <- COLOR_THEME_MAPPING[[key]]
    if (is.null(input_ids)) next
    val <- vals[[key]]
    if (is.null(val) || !nzchar(as.character(val)[1])) next
    for (s in targets) {
      for (inp in input_ids) {
        js <- c(js, sprintf(
          "Shiny.setInputValue('%s','%s');", s$ns(inp), val
        ))
        colourpicker::updateColourInput(s, inp, value = val)
      }
    }
  }
  if (length(js)) {
    shinyjs::runjs(paste(js, collapse = ""))
  }
  invisible(NULL)
}
