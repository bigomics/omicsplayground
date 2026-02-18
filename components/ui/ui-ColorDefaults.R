##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## ---------------------------------------------------------------
## Color Theme Defaults
## ---------------------------------------------------------------

COLOR_THEME_DEFAULTS <- list(
  primary    = "#f23451",   # Up / High
  secondary  = "#3181de",   # Down / Low / scatter_color
  neutral    = "#eeeeee",   # Mid (heatmap midpoint)
  ns_color   = "#eeeeee",   # Not significant dots
  bar_color  = "#A6CEE3",   # Bar color
  accent     = "#e3a45a",   # Significant in one (color_one)
  success    = "#5B9B5B",   # Significant in both (color_both)
  line       = "#00EE00",   # Enrichment line
  palette    = "default",
  palette_c1 = "#3181de",   # Custom gradient start
  palette_c2 = "#eeeeee",   # Custom gradient middle
  palette_c3 = "#f23451"    # Custom gradient end
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
## Singleton reactive store (module-level environment)
## ---------------------------------------------------------------

.color_theme_env <- new.env(parent = emptyenv())
.color_theme_env$theme <- NULL

#' Initialise the colour theme reactive store (call once in server).
#' @return invisible reactiveValues
init_color_theme <- function() {
  if (is.null(.color_theme_env$theme)) {
    .color_theme_env$theme <- do.call(shiny::reactiveValues, COLOR_THEME_DEFAULTS)
  }
  invisible(.color_theme_env$theme)
}

#' Return the colour theme reactive store.
#' @return reactiveValues
get_color_theme <- function() {
  if (is.null(.color_theme_env$theme)) {
    init_color_theme()
  }
  .color_theme_env$theme
}
