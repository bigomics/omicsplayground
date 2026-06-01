# copilot_pgx_normalize.R — Centralised PGX normaliser
#
# Single funnel that every code path uses to convert whatever pgx-like
# value it has (NULL, reactiveValues, plain list, classed pgx list) into
# a value that satisfies omicspgx::.pgx_check() — i.e. a plain list with
# class c("pgx", "list") — or NULL when no data is available.
#
# Why this exists
# ---------------
# The chat board receives PGX from three sources:
#   1. Global `pgx` reactiveValues observer (CopilotBoardServer.R)
#   2. `pgx_loaded_event` from the manage_pgx tool (omicsagentovi)
#   3. The currently-loaded agent's `@context@pgx` slot
#   4. The restore controller, reading the global pgx at restore time
#
# Historically each site re-implemented the snapshot+class stamp by hand,
# and at least one (the restore controller) forgot — leaking the bare
# reactiveValues into agent@context@pgx and tripping .pgx_check() with
# "got 'reactivevalues'" on the next omicspgx call (e.g. after clicking
# "+ New chat" the freshly rebuilt agent inherited the dirty pgx).
#
# Funnelling every site through copilot_normalize_pgx() guarantees:
#   - reactiveValues are isolated + snapshot-listed exactly once
#   - the class stamp is always applied
#   - one place to add tracing or extra coercions
#' Normalise a pgx-shaped value into a classed pgx list
#'
#' Accepts `NULL`, a `shiny::reactiveValues`, a plain `list`, or an
#' already classed `pgx` list. Returns a list with class
#' `c("pgx", "list")` or `NULL` when the input has no usable payload.
#'
#' @param pgx Any pgx-shaped value (or `NULL`).
#' @param source Character tag identifying the caller. Reserved for future
#'   tracing; ignored by the current implementation.
#' @return A classed pgx list, or `NULL`.
#' @export
copilot_normalize_pgx <- function(pgx, source = "unknown") {
  if (is.null(pgx)) return(NULL)

  # reactiveValues -> plain list. Must isolate; reactiveValuesToList()
  # touches the dependency graph and will explode outside a reactive
  # consumer (e.g. session$onFlushed callbacks, observer bodies that
  # already isolated upstream).
  #
  # Log every actual coercion: this is supposed to be a fallback for
  # source-1/2 boundaries that already snapshot, so if the call site is
  # anything else (e.g. `run_controller/.current_pgx/agent`) we have a
  # regression somewhere upstream that we want to see in the logs.
  if (inherits(pgx, "reactivevalues")) {
    if (exists("log_info", mode = "function")) {
      log_info("copilot.pgx.normalize_coerce",
               source = source %||% "unknown")
    }
    pgx <- tryCatch(
      shiny::isolate(shiny::reactiveValuesToList(pgx)),
      error = function(e) NULL
    )
    if (is.null(pgx)) return(NULL)
  }

  # Sanity: must look like a populated pgx (has $X and $name). Otherwise
  # treat as "no dataset yet" and return NULL so downstream first-load
  # guards (apply_dataset, .current_pgx) keep behaving the same way.
  if (is.null(pgx[["X"]]) || is.null(pgx[["name"]])) return(NULL)

  # Stamp the class. unique() preserves an existing "pgx" tag from inputs
  # already produced by playbase::pgx.initialize().
  if (!inherits(pgx, "pgx")) {
    class(pgx) <- unique(c("pgx", class(pgx)))
  }

  pgx
}
