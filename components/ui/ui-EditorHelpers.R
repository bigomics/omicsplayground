##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## ---------------------------------------------------------------
## Editor Helper Functions
## ---------------------------------------------------------------
## Centralised helpers that eliminate repeated boilerplate across
## 60+ plot editor modules.  Sourced via 00SourceAll.R.
## ---------------------------------------------------------------


## ---------------------------------------------------------------
## 1. Color helpers
## ---------------------------------------------------------------

#' Resolve an editor colour input with fallback.
#'
#' @param input   Shiny input object.
#' @param id      Character input ID (e.g. "color_up").
#' @param default Fallback: a theme key like "primary" (looked up via
#'   \code{get_color_theme()}) or a hex colour like "#f23451" (used as-is).
#' @return A single colour string.
get_editor_color <- function(input, id, default) {
  val <- input[[id]]
  if (!is.null(val)) return(val)
  if (grepl("^#", default)) default else get_color_theme()[[default]]
}

#' Build the named colour vector used by volcano / MA-plot renderers.
#'
#' @param input  Shiny input list.
#' @param notsig Colour for non-significant points (default \code{"#707070AA"}).
#' @param notsel Colour for non-selected points (default \code{"#cccccc88"}).
#'   Set to \code{NULL} to omit.
#' @return A named character vector.
extract_volcano_colors <- function(input,
                                   notsig = "#707070AA",
                                   notsel = "#cccccc88") {
  col_up   <- get_editor_color(input, "color_up",   "primary")
  col_down <- get_editor_color(input, "color_down", "secondary")
  colors <- c(up = col_up, notsig = notsig, down = col_down)
  if (!is.null(notsel)) colors <- c(colors, notsel = notsel)
  colors
}

#' Build the three-element heatmap colour vector (low, mid, high).
#'
#' Returns \code{c(col_down, mid, col_up)} matching the positional
#' convention used by \code{pgx.plotActivation()} and friends.
#'
#' @param input Shiny input list.
#' @param mid   Default midpoint colour (default \code{"grey90"}).
#' @return Character vector of length 3.
extract_heatmap_colors <- function(input, mid = "grey90") {
  col_up   <- get_editor_color(input, "color_up",   "primary")
  col_down <- get_editor_color(input, "color_down", "secondary")
  c(col_down, mid, col_up)
}

#' Build heatmap colours from separate low/mid/high editor inputs.
#'
#' Used by editors that expose \code{color_low}, \code{color_mid},
#' \code{color_high} (e.g. split-heatmap).
#'
#' @param input Shiny input list.
#' @return Character vector of length 3.
extract_heatmap_lmh_colors <- function(input) {
  c(
    get_editor_color(input, "color_low",  "secondary"),
    get_editor_color(input, "color_mid",  "neutral"),
    get_editor_color(input, "color_high", "primary")
  )
}


## ---------------------------------------------------------------
## 2. Custom label resolution
## ---------------------------------------------------------------

#' Resolve custom label features from the plot editor.
#'
#' Replaces the repeated pattern of checking \code{input$custom_labels},
#' calling \code{parse_label_features()}, and optionally running
#' \code{map_probes()}.
#'
#' @param input         Shiny input object.
#' @param feature_names Character vector of valid feature names to match against.
#' @param defaults      Value to return when custom labelling is off (can be NULL).
#' @param pgx           Optional pgx object; when provided, results are passed
#'   through \code{playbase::map_probes(pgx$genes, ...)}.
#' @return Character vector of resolved labels, or \code{defaults}.
get_custom_labels <- function(input, feature_names, defaults = NULL, pgx = NULL) {
  if (!isTRUE(input$custom_labels)) return(defaults)

  label_text <- input$label_features
  if (is.null(label_text) || !nzchar(trimws(label_text))) return(defaults)

  resolved <- parse_label_features(label_text, feature_names)
  if (is.null(resolved) || length(resolved) == 0) return(defaults)

  if (!is.null(pgx)) {
    resolved <- playbase::map_probes(pgx$genes, resolved)
  }

  resolved
}


## ---------------------------------------------------------------
## 3. Margin & aspect ratio (ggplot2)
## ---------------------------------------------------------------

#' Apply editor margin and aspect-ratio theme to a ggplot.
#'
#' Reads \code{margin_checkbox}, \code{margin_top/right/bottom/left},
#' \code{aspect_ratio_checkbox}, and \code{aspect_ratio} from Shiny
#' \code{input} and applies them as \code{ggplot2::theme()} layers.
#'
#' @param p     A ggplot2 plot object.
#' @param input Shiny input list.
#' @param margin_default       Default margin in pt (default 10).
#' @param aspect_ratio_default Default aspect ratio (default 0.5).
#' @return The (possibly modified) ggplot2 plot object.
apply_editor_theme <- function(p, input,
                               margin_default = 10,
                               aspect_ratio_default = 0.5) {
  if (isTRUE(input$margin_checkbox)) {
    mt <- ifelse(is.na(input$margin_top),    margin_default, input$margin_top)
    mr <- ifelse(is.na(input$margin_right),  margin_default, input$margin_right)
    mb <- ifelse(is.na(input$margin_bottom), margin_default, input$margin_bottom)
    ml <- ifelse(is.na(input$margin_left),   margin_default, input$margin_left)
    p <- p + ggplot2::theme(
      plot.margin = ggplot2::margin(t = mt, r = mr, b = mb, l = ml, unit = "pt")
    )
  }

  if (isTRUE(input$aspect_ratio_checkbox)) {
    ar <- if (is.na(input$aspect_ratio)) aspect_ratio_default else input$aspect_ratio
    p <- p + ggplot2::theme(aspect.ratio = ar)
  }

  p
}


## ---------------------------------------------------------------
## 4. ggprism settings extraction
## ---------------------------------------------------------------

#' Extract ggprism settings from Shiny input.
#'
#' Returns a named list that can be spliced into \code{ggVolcano()} or
#' similar calls via \code{do.call()}.
#'
#' @param input Shiny input object.
#' @return Named list of ggprism parameters.
extract_ggprism_params <- function(input) {
  list(
    use_ggprism           = isTRUE(input$use_ggprism),
    ggprism_palette       = if (is.null(input$ggprism_palette)) "black_and_white" else input$ggprism_palette,
    ggprism_colors        = isTRUE(input$ggprism_colors),
    ggprism_border        = isTRUE(input$ggprism_border),
    ggprism_axis_guide    = if (is.null(input$ggprism_axis_guide)) "default" else input$ggprism_axis_guide,
    ggprism_show_legend   = isTRUE(input$ggprism_show_legend),
    ggprism_legend_x      = if (is.null(input$ggprism_legend_x) || is.na(input$ggprism_legend_x)) 0.95 else input$ggprism_legend_x,
    ggprism_legend_y      = if (is.null(input$ggprism_legend_y) || is.na(input$ggprism_legend_y)) 0.95 else input$ggprism_legend_y,
    ggprism_legend_border = isTRUE(input$ggprism_legend_border)
  )
}


## ---------------------------------------------------------------
## 5. Label rendering settings extraction
## ---------------------------------------------------------------

#' Extract label rendering settings from Shiny input.
#'
#' Reads optional plot-editor knobs (label size, marker size, etc.)
#' with sensible defaults.
#'
#' @param input    Shiny input object.
#' @param defaults Optional named list to override built-in defaults
#'   (e.g. \code{list(label_size = 4, marker_size = 1)}).
#' @return Named list of label settings.
extract_label_settings <- function(input, defaults = list()) {
  dfl <- list(
    label_size         = 5,
    marker_size        = 1.2,
    axis_text_size     = 14,
    box_padding        = 0.1,
    min_segment_length = 0,
    label_box          = TRUE,
    segment_linetype   = 1L,
    hyperbola_k        = 1
  )
  dfl[names(defaults)] <- defaults

  ## Safely read one input value, guarding against NULL and optionally NA.
  safe <- function(key, default, null_na = FALSE, cast = identity) {
    val <- input[[key]]
    missing <- is.null(val) || (null_na && isTRUE(is.na(val)))
    if (missing) default else cast(val)
  }

  list(
    label_size         = safe("label_size",         dfl$label_size),
    marker_size        = safe("marker_size",        dfl$marker_size),
    axis_text_size     = safe("axis_text_size",     dfl$axis_text_size),
    box_padding        = safe("box_padding",        dfl$box_padding,        null_na = TRUE),
    min_segment_length = safe("min_segment_length", dfl$min_segment_length, null_na = TRUE),
    label_box          = safe("label_box",          dfl$label_box),
    segment_linetype   = safe("segment_linetype",   dfl$segment_linetype,   cast = as.integer),
    hyperbola_k        = safe("hyperbola_k",        dfl$hyperbola_k)
  )
}


## ---------------------------------------------------------------
## 6. Custom palette UI helpers
## ---------------------------------------------------------------

#' Generate colour-picker inputs for a custom palette.
#'
#' Creates one \code{colourpicker::colourInput} per label, using
#' \code{custom_color_1}, \code{custom_color_2}, ... as input IDs.
#' Call this inside \code{output$custom_palette_ui <- renderUI(...)}.
#'
#' @param labels       Character vector of group/series names to label each picker.
#' @param ns_func      Namespace function (\code{ns} or \code{session$ns}).
#' @param default_colors Optional character vector of default hex colours.
#'   Recycled if shorter than \code{labels}. When \code{NULL}, uses
#'   \code{omics_pal_d("muted_light")}.
#' @return A \code{shiny::tagList} of colourInput widgets.
custom_palette_pickers <- function(labels, ns_func,
                                   default_colors = NULL) {
  n <- length(labels)
  if (is.null(default_colors)) {
    default_colors <- omics_pal_d("muted_light")(n)
  }
  ## Recycle colours if fewer than labels
  default_colors <- rep_len(default_colors, n)

  pickers <- lapply(seq_along(labels), function(i) {
    colourpicker::colourInput(
      ns_func(paste0("custom_color_", i)),
      label = labels[i],
      value = default_colors[i]
    )
  })
  shiny::tagList(pickers)
}

#' Read back custom palette colours from Shiny input.
#'
#' Reads \code{input$custom_color_1} .. \code{input$custom_color_n},
#' falling back to \code{fallback_colors} for any that are NULL
#' (not yet rendered).
#'
#' @param input           Shiny input object.
#' @param n               Number of colours to read.
#' @param fallback_colors Default colours when an input is NULL.
#'   Recycled if shorter than \code{n}. When \code{NULL}, uses
#'   \code{omics_pal_d("muted_light")}.
#' @return Character vector of length \code{n}.
get_custom_palette_colors <- function(input, n, fallback_colors = NULL) {
  if (is.null(fallback_colors)) {
    fallback_colors <- omics_pal_d("muted_light")(n)
  }
  fallback_colors <- rep_len(fallback_colors, n)
  sapply(seq_len(n), function(j) {
    val <- input[[paste0("custom_color_", j)]]
    if (is.null(val)) fallback_colors[j] else val
  })
}


## ---------------------------------------------------------------
## 7. Sortable rank list UI
## ---------------------------------------------------------------

#' Create a sortable bucket list for drag-and-drop ordering.
#'
#' Wraps \code{sortable::bucket_list} + \code{sortable::add_rank_list}
#' with the standard editor styling.  Call inside
#' \code{output$rank_list <- renderUI(...)}.
#'
#' @param labels   Character vector of items to display / reorder.
#' @param ns_func  Namespace function (\code{ns} or \code{session$ns}).
#' @param input_id Character input ID for the rank list result
#'   (default \code{"rank_list_basic"}).
#' @return A \code{sortable::bucket_list} tag.
rank_list_ui <- function(labels, ns_func, input_id = "rank_list_basic") {
  sortable::bucket_list(
    header = NULL,
    class = "default-sortable custom-sortable",
    sortable::add_rank_list(
      input_id = ns_func(input_id),
      text = NULL,
      labels = labels
    )
  )
}
