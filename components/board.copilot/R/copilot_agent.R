#' Copilot Agent Backend
#'
#' Wraps omicsagentovi::Agent for use in the copilot board.
#' Returns a list with $prompt(text), $agent, $reset(), $set_pgx(pgx)

copilot_as_pgx <- function(pgx) {
  if (is.null(pgx)) {
    return(NULL)
  }

  if (!inherits(pgx, "pgx")) {
    class(pgx) <- unique(c("pgx", class(pgx)))
  }

  pgx
}

copilot_snapshot_pgx <- function(pgx) {
  if (is.null(pgx) || is.null(pgx$X)) {
    return(NULL)
  }

  snapshot <- shiny::reactiveValuesToList(pgx)
  copilot_as_pgx(snapshot)
}

copilot_create_context <- function(pgx = NULL, plot_callback = NULL, pgx_dir = NULL) {
  ctx <- omicsagentovi::ovi_context(copilot_as_pgx(pgx))

  if (is.function(plot_callback)) {
    ctx@state$plot_callback <- plot_callback
  }
  if (!is.null(pgx_dir)) {
    ctx@state$data_dir <- pgx_dir
  }

  ctx
}

copilot_create_agent <- function(pgx = NULL, plot_callback = NULL,
                                  pgx_dir = NULL, tier = "copilot-default") {
  info("[copilot_create_agent] creating ovi_context; tier=", tier)
  ctx <- copilot_create_context(
    pgx = pgx,
    plot_callback = plot_callback,
    pgx_dir = pgx_dir
  )

  agent <- omicsagentovi::Agent(tier = tier, context = ctx)
  model_name <- tryCatch(agent@model, error = function(e) "unknown")
  info("[copilot_create_agent] Agent created — tier=", tier, " model=", model_name)

  list(
    prompt = function(text) {
      omicsagentovi::agent_prompt(agent, text)
    },
    agent = agent,
    tier = tier,
    set_pgx = function(new_pgx) {
      omicsagentovi::ovi_set_pgx(agent, copilot_as_pgx(new_pgx))
    },
    reset = function(new_pgx = NULL, new_plot_callback = NULL, new_tier = NULL) {
      info("[copilot_create_agent] reset called")
      pgx_to_use      <- if (is.null(new_pgx)) pgx else new_pgx
      callback_to_use <- if (is.function(new_plot_callback)) new_plot_callback else plot_callback
      tier_to_use     <- if (!is.null(new_tier)) new_tier else tier
      ctx   <<- copilot_create_context(
        pgx = pgx_to_use,
        plot_callback = callback_to_use,
        pgx_dir = pgx_dir
      )
      agent <<- omicsagentovi::Agent(tier = tier_to_use, context = ctx)
      tier  <<- tier_to_use
      model_name <- tryCatch(agent@model, error = function(e) "unknown")
      info("[copilot_create_agent] Agent recreated — tier=", tier_to_use, " model=", model_name)
    }
  )
}
