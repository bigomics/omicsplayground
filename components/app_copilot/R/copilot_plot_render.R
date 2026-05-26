# copilot_plot_render.R — Recipe renderer for Copilot evidence plots
#
# Pure helpers: no reactive dependencies, no Shiny state.
# This IS the board-side recipe renderer. The package (omicsagentovi) delivers
# recipes (plot_type + args); this module executes them.
#
# Called from the plot_callback closure in copilot_bindings.R, which receives
# (pgx, plot_type, args, artifact) from the package and calls:
#   1. copilot_build_plot(pgx, plot_type, args)   -> R plot object
#   2. copilot_detect_plot_kind(plot_obj)          -> "ggplot"|"plotly"|"iheatmapr"|NULL
#   3. copilot_prerender_ggplot/copilot_prerender_plotly as appropriate
#   4. evidence_api$append_artifact(record)

# ---- Null-coalescing operator (local — avoids dependency on run_controller.R) ----
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

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
#' Dispatches to the appropriate omicspgxmcp / playbase renderer based
#' on \code{plot_type}. The PGX must already be class-tagged (plain list
#' with class "pgx") — coercion via `copilot_as_pgx` has been removed
#' because the agent context guarantees this at the boundary.
#'
#' @param pgx A PGX-shaped list (class already set by apply_dataset boundary).
#' @param plot_type Character scalar — one of "pca", "tsne", "volcano",
#'   "heatmap", "ma", "barplot_de", "enrichment_dotplot".
#' @param args Named list with optional \code{contrast}, \code{features},
#'   and \code{collection} entries.
#' @return A plot object (ggplot, plotly, or iheatmapr).
copilot_build_plot <- function(pgx, plot_type, args) {
  # NOTE: copilot_as_pgx() call removed — PGX is already class-tagged at
  # the agent context boundary (copilot_server.R global observer + run
  # controller .current_pgx). Calling it here would double-coerce.
  plot_type  <- tolower(trimws(plot_type))
  contrast   <- if (!is.null(args$contrast))   args$contrast   else NULL
  features   <- copilot_parse_features(if (!is.null(args$features)) args$features else NULL)
  collection <- if (!is.null(args$collection)) args$collection else NULL

  switch(plot_type,
    pca = {
      params    <- omicspgxmcp:::new_plot_params("scatter", params = list(method = "pca"))
      extracted <- omicspgxmcp:::.extract_scatter_data(pgx, params)
      do.call(omicsplots::pgx.plot_scatter, extracted$renderer_args)
    },
    tsne = {
      params    <- omicspgxmcp:::new_plot_params("scatter", params = list(method = "tsne"))
      extracted <- omicspgxmcp:::.extract_scatter_data(pgx, params)
      do.call(omicsplots::pgx.plot_scatter, extracted$renderer_args)
    },
    volcano = {
      params    <- omicspgxmcp:::new_plot_params("volcano", contrast = contrast, params = list())
      extracted <- omicspgxmcp:::.extract_volcano_data(pgx, params)
      do.call(
        omicsplots::pgx.plot_volcano,
        c(extracted$renderer_args, list(show_sample_badge = FALSE))
      )
    },
    heatmap = {
      params    <- omicspgxmcp:::new_plot_params(
        "heatmap",
        contrast = contrast,
        params   = list(genes = features)
      )
      extracted <- omicspgxmcp:::.extract_heatmap_data(pgx, params)
      do.call(omicsplots::pgx.plot_heatmap, extracted$renderer_args)
    },
    ma = {
      playbase::pgx.plotMA(pgx, contrast = contrast, plotlib = "ggplot")
    },
    barplot_de = {
      params    <- omicspgxmcp:::new_plot_params(
        "barplot",
        contrast = contrast,
        params   = list(what = "de")
      )
      extracted <- omicspgxmcp:::.extract_barplot_data(pgx, params)
      do.call(omicsplots::pgx.plot_barplot, extracted$renderer_args)
    },
    enrichment_dotplot = {
      playbase::pgx.plotEnrichmentDotPlot(pgx, contrast = contrast, filter = collection)
    },
    stop(sprintf("Unsupported copilot plot type: %s", plot_type), call. = FALSE)
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
