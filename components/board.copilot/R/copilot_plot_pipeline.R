#' Copilot Plot Pipeline
#'
#' Pure helper functions for building, pre-rendering, and managing
#' copilot evidence plots.  No reactive dependencies — safe to call
#' from any context.

#' Detect the plot kind from an R plot object
#'
#' @param plot_obj Any R object returned by an agent visualization tool.
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
#' on \code{plot_type}.
#'
#' @param pgx A PGX-shaped list (will be coerced via \code{copilot_as_pgx}).
#' @param plot_type Character scalar — one of "pca", "tsne", "volcano",
#'   "heatmap", "ma", "barplot_de", "enrichment_dotplot".
#' @param args Named list with optional \code{contrast}, \code{features},
#'   and \code{collection} entries.
#' @return A plot object (ggplot, plotly, or iheatmapr).
copilot_build_plot <- function(pgx, plot_type, args) {
  pgx <- copilot_as_pgx(pgx)
  plot_type <- tolower(trimws(plot_type))
  contrast <- if (!is.null(args$contrast)) args$contrast else NULL
  features <- copilot_parse_features(if (!is.null(args$features)) args$features else NULL)
  collection <- if (!is.null(args$collection)) args$collection else NULL

  switch(plot_type,
    pca = {
      params <- omicspgxmcp:::new_plot_params("scatter", params = list(method = "pca"))
      extracted <- omicspgxmcp:::.extract_scatter_data(pgx, params)
      do.call(omicsplots::pgx.plot_scatter, extracted$renderer_args)
    },
    tsne = {
      params <- omicspgxmcp:::new_plot_params("scatter", params = list(method = "tsne"))
      extracted <- omicspgxmcp:::.extract_scatter_data(pgx, params)
      do.call(omicsplots::pgx.plot_scatter, extracted$renderer_args)
    },
    volcano = {
      params <- omicspgxmcp:::new_plot_params("volcano", contrast = contrast, params = list())
      extracted <- omicspgxmcp:::.extract_volcano_data(pgx, params)
      do.call(
        omicsplots::pgx.plot_volcano,
        c(extracted$renderer_args, list(show_sample_badge = FALSE))
      )
    },
    heatmap = {
      params <- omicspgxmcp:::new_plot_params(
        "heatmap",
        contrast = contrast,
        params = list(genes = features)
      )
      extracted <- omicspgxmcp:::.extract_heatmap_data(pgx, params)
      do.call(omicsplots::pgx.plot_heatmap, extracted$renderer_args)
    },
    ma = {
      playbase::pgx.plotMA(
        pgx,
        contrast = contrast,
        plotlib = "ggplot"
      )
    },
    barplot_de = {
      params <- omicspgxmcp:::new_plot_params(
        "barplot",
        contrast = contrast,
        params = list(what = "de")
      )
      extracted <- omicspgxmcp:::.extract_barplot_data(pgx, params)
      do.call(omicsplots::pgx.plot_barplot, extracted$renderer_args)
    },
    enrichment_dotplot = {
      playbase::pgx.plotEnrichmentDotPlot(
        pgx,
        contrast = contrast,
        filter = collection
      )
    },
    stop(sprintf("Unsupported copilot plot type: %s", plot_type), call. = FALSE)
  )
}

#' FIFO-prune persisted chat sessions to a maximum count
#'
#' Deletes the oldest sessions (by \code{updated_at}) when the total
#' exceeds \code{max_sessions}.
#'
#' @param store A \code{SessionStore} object.
#' @param max_sessions Integer maximum number of sessions to keep.
#' @return Invisible NULL.
copilot_prune_sessions <- function(store, max_sessions = 100L) {
  session_dir <- store@state$session_dir
  ids <- omicsagentovi::ovi_sessions(session_dir = session_dir)
  if (length(ids) <= max_sessions) return(invisible(NULL))
  metas <- lapply(ids, function(i) {
    tryCatch(
      omicsagentovi::ovi_session_meta(i, session_dir = session_dir),
      error = function(e) NULL
    )
  })
  updated <- vapply(metas, function(m) m$updated_at %||% "", character(1L))
  ord <- order(updated)  # oldest first
  n_drop <- length(ids) - max_sessions
  for (id in ids[ord[seq_len(n_drop)]]) {
    try(unlink(file.path(session_dir, id), recursive = TRUE), silent = TRUE)
  }
  invisible(NULL)
}

#' Replay persisted user/assistant text turns into a shinychat instance
#'
#' Tool turns are intentionally skipped per the no-tool-rendering policy.
#'
#' @param chat_id The shinychat widget id to append messages to.
#' @param turns List of ellmer Turn objects from a restored session.
#' @return Integer count of user turns replayed.
copilot_replay_turns <- function(chat_id, turns) {
  n_user <- 0L
  for (turn in turns) {
    role <- tryCatch(turn@role, error = function(e) NULL)
    if (is.null(role) || !(role %in% c("user", "assistant"))) next
    contents <- tryCatch(turn@contents, error = function(e) list())
    has_tool <- any(vapply(contents, function(c) {
      S7::S7_inherits(c, ellmer::ContentToolRequest) ||
      S7::S7_inherits(c, ellmer::ContentToolResult)
    }, logical(1)))
    if (has_tool) next
    text_parts <- vapply(contents, function(c) {
      if (S7::S7_inherits(c, ellmer::ContentText)) c@text %||% "" else ""
    }, character(1))
    text <- paste(text_parts[nzchar(text_parts)], collapse = "\n")
    if (!nzchar(text)) next
    shinychat::chat_append_message(
      chat_id, list(role = role, content = text), chunk = FALSE
    )
    if (role == "user") n_user <- n_user + 1L
  }
  n_user
}

#' Human-readable label for a copilot tier identifier
#'
#' @param tier Character scalar tier id (e.g. "copilot-default").
#' @return Character scalar display label.
copilot_tier_label <- function(tier) {
  switch(tier,
    `copilot-default` = "Balanced (Default)",
    `copilot-fast`    = "Fast",
    `copilot-deep`    = "Deep Think",
    tier
  )
}

#' Pre-render a ggplot to a temporary PNG file
#'
#' Renders a ggplot via \code{ggplot2::ggsave()} and returns the file
#' path.  The cost of the graphics device and PNG encoding happens here,
#' outside \code{flushReact()}.
#'
#' @param plot_obj A ggplot object.
#' @param width Plot width in inches (default 8).
#' @param height Plot height in inches (default 6).
#' @param dpi Resolution in dots per inch (default 96).
#' @return Character path to the temporary PNG file.
copilot_prerender_ggplot <- function(plot_obj, width = 8, height = 6, dpi = 96) {
  path <- tempfile(fileext = ".png")
  ggplot2::ggsave(path, plot = plot_obj, width = width, height = height, dpi = dpi)
  dbg("[CopilotPrerender] ggplot -> PNG path=", path, " size=", file.size(path))
  path
}

#' Pre-build a plotly object for zero-overhead rendering
#'
#' Calls \code{plotly::plotly_build()} and strips unnecessary modebar
#' buttons so that \code{renderPlotly} can relay the object with no
#' JSON serialization or layout computation cost.
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
  dbg("[CopilotPrerender] plotly_build done")
  built
}

#' Clean up pre-rendered temp files
#'
#' Removes temporary PNG files produced by \code{copilot_prerender_ggplot()}.
#' Silently ignores missing files and unlink errors.
#'
#' @param paths Character vector of file paths (NULLs are tolerated).
#' @return Invisible NULL.
copilot_prerender_cleanup <- function(paths) {
  for (p in paths) {
    if (!is.null(p) && file.exists(p)) {
      tryCatch(unlink(p), error = function(e) NULL)
    }
  }
  invisible(NULL)
}
