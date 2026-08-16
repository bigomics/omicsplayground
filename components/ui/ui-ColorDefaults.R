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

  theme <- get_color_theme()

  ## For each theme key, observe changes and push to editor colour inputs.
  ## We use Shiny.setInputValue via JS to directly set the server-side
  ## input value.  This is more reliable than updateColourInput because
  ## the colourpicker widget may not be initialised when its container
  ## (the editor modal) is hidden.  We also call updateColourInput so
  ## that when the user does open the editor modal, the widget shows the
  ## correct swatch colour.
  ## ignoreInit = FALSE so lazily-loaded modules pick up values that
  ## were changed before the module was navigated to.
  lapply(names(COLOR_THEME_MAPPING), function(key) {
    input_ids <- COLOR_THEME_MAPPING[[key]]
    shiny::observeEvent(theme[[key]],
      {
        val <- theme[[key]]
        for (inp in input_ids) {
          full_id <- parent_session$ns(inp)
          shinyjs::runjs(sprintf(
            "Shiny.setInputValue('%s', '%s')",
            full_id, val
          ))
          colourpicker::updateColourInput(parent_session, inp, value = val)
        }
      },
      ignoreInit = FALSE
    )
  })

  ## Palette observer
  shiny::observeEvent(theme$palette,
    {
      shiny::updateSelectInput(parent_session, "palette", selected = theme$palette)
    },
    ignoreInit = FALSE
  )
}
