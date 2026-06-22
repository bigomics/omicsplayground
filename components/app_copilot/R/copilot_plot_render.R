# copilot_plot_render.R — Recipe renderer for Copilot evidence plots
#
# Pure helpers: no reactive dependencies, no Shiny state.
# This is the board-side plot result adapter. The package (omicsagentovi)
# normally delivers a prebuilt plot_result; recipe fallback is normalized into
# the same omicspgxmcp shared plot contract.
#
# Called from the plot_callback closure in copilot_bindings.R, which receives
# (pgx, plot_type, args, artifact) from the package and calls:
#   1. plot_result$plot or copilot_build_plot(pgx, plot_type, args)
#   2. copilot_detect_plot_kind(plot_obj)          -> "ggplot"|"plotly"|"iheatmapr"|NULL
#   3. copilot_prerender_ggplot/copilot_prerender_plotly as appropriate
#   4. evidence_api$append_artifact(record)

#' Detect the plot kind from an R plot object
#'
#' @param plot_obj Any R object returned by an agent visualisation tool.
#' @return Character scalar: `"ggplot"`, `"plotly"`, `"iheatmapr"`, or `NULL`.
copilot_detect_plot_kind <- function(plot_obj) {
  if (inherits(plot_obj, "plotly"))    return("plotly")
  if (inherits(plot_obj, "iheatmapr")) return("iheatmapr")
  if (inherits(plot_obj, "ggplot"))    return("ggplot")
  NULL
}

#' Parse a comma-separated feature string into a character vector
#'
#' @param features Comma-separated string of feature names (or NULL).
#' @return Character vector of trimmed, non-empty feature names, or NULL.
copilot_parse_features <- function(features) {
  if (is.null(features) || !nzchar(trimws(features))) {
    return(NULL)
  }
  out <- trimws(unlist(strsplit(features, ",")))
  out[nzchar(out)]
}

#' Build a copilot evidence plot from a PGX object
#'
#' Normalizes older webapp recipe names into the shared omicspgxmcp plot
#' contract and delegates all plot construction to `pgx_build_omics_plot()`.
#'
#' @param pgx A PGX-shaped list (class already set by apply_dataset boundary).
#' @param plot_type Character scalar plot type or compatibility alias.
#' @param args Named list with `params`, `contrast`, or older recipe fields.
#' @return A plot object (ggplot, plotly, or iheatmapr).
copilot_build_plot <- function(pgx, plot_type, args) {
  request <- copilot_normalize_plot_request(plot_type, args)
  omicspgxmcp::pgx_build_omics_plot(
    pgx = pgx,
    plot_type = request$plot_type,
    params = request$params,
    target = request$target,
    include_renderer_args = FALSE
  )$plot
}

copilot_normalize_plot_request <- function(plot_type, args) {
  raw_type <- tolower(trimws(plot_type))
  params <- args$params %||% list()
  if (!is.list(params)) params <- as.list(params)

  if (!is.null(args$contrast) && is.null(params$contrast)) {
    params$contrast <- args$contrast
  }
  if (!is.null(args$title) && is.null(params$title)) {
    params$title <- args$title
  }

  if (identical(raw_type, "pca")) {
    plot_type <- "scatter"
    params$method <- "pca"
  } else if (identical(raw_type, "tsne")) {
    plot_type <- "scatter"
    params$method <- "tsne"
  } else if (identical(raw_type, "barplot_de")) {
    plot_type <- "barplot"
    params$what <- "de"
  } else if (identical(raw_type, "enrichment_dotplot")) {
    plot_type <- "dotplot"
    params$what <- "enrichment"
    if (!is.null(params$collection) && is.null(params$query)) {
      params$query <- params$collection
    }
    params$collection <- NULL
    if (!is.null(args$collection) && is.null(params$query)) {
      params$query <- args$collection
    }
  } else if (identical(raw_type, "corrplot")) {
    plot_type <- "corrmat"
  } else if (identical(raw_type, "maplot")) {
    plot_type <- "ma"
  } else {
    plot_type <- raw_type
  }

  if (!is.null(args$features) && is.null(params$genes)) {
    params$genes <- copilot_parse_features(args$features)
  }

  list(
    plot_type = plot_type,
    params = params,
    target = args$target %||% "user"
  )
}

#' Pre-render a ggplot to a temporary PNG file
#'
#' Renders the plot via \code{ggplot2::ggsave()} and returns the file path.
#' The cost of the graphics device and PNG encoding happens here, outside
#' the reactive flush cycle.
#'
#' @param plot_obj A ggplot object.
#' @param width Plot width in inches (default 8).
#' @param height Plot height in inches (default 6).
#' @param dpi Resolution in dots per inch (default 96).
#' @return Character path to the temporary PNG file.
copilot_prerender_ggplot <- function(plot_obj, width = 8, height = 6, dpi = 96) {
  path <- tempfile(fileext = ".png")
  ggplot2::ggsave(path, plot = plot_obj, width = width, height = height, dpi = dpi)
  log_info("copilot.prerender.ggplot", path = path, size = file.size(path))
  path
}

#' Pre-build a plotly object for zero-overhead rendering
#'
#' Calls \code{plotly::plotly_build()} and strips unnecessary modebar buttons
#' so that \code{renderPlotly} can relay the object with no JSON serialisation
#' or layout computation cost.
#'
#' @param plot_obj A plotly object.
#' @return A pre-built plotly object ready for \code{renderPlotly}.
copilot_prerender_plotly <- function(plot_obj) {
  built <- plotly::plotly_build(plot_obj)
  built <- plotly::config(
    built,
    modeBarButtonsToRemove = list(
      "sendDataToCloud", "editInChartStudio", "lasso2d", "select2d",
      "autoScale2d", "hoverClosestCartesian", "hoverCompareCartesian",
      "toggleSpikelines"
    ),
    modeBarButtonsToAdd = list(),
    displaylogo = FALSE
  )
  # Register plotly_click so the global PlotModule click listener
  # (event_data("plotly_click") in ui-PlotModule.R) doesn't warn about an
  # unregistered source "A" when this plot is the one being clicked.
  # CAVEAT: that PlotModule observer also reads click_data$key and writes to
  # parent_session$input$label_features. Because copilot plots inherit the
  # default plotly source "A", clicks on agent traces carrying a `key` field
  # will leak into a board's "Label features" textarea. Fix when it bites:
  # set a non-default source on copilot plots, or scope the PlotModule
  # observer to a specific source.
  built <- plotly::event_register(built, "plotly_click")
  log_info("copilot.prerender.plotly")
  built
}

#' Clean up pre-rendered temp PNG files
#'
#' Removes temporary PNG files produced by \code{copilot_prerender_ggplot()}.
#' Silently ignores missing files and unlink errors. NULL/NA entries tolerated.
#'
#' @param paths Character vector of file paths (NULLs and NAs tolerated).
#' @return Invisible NULL.
copilot_prerender_cleanup <- function(paths) {
  for (p in paths) {
    if (!is.null(p) && !is.na(p) && file.exists(p)) {
      tryCatch(unlink(p), error = function(e) NULL)
    }
  }
  invisible(NULL)
}
