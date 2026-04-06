#' Copilot Board Server
#'
#' Agent lifecycle, chat wiring, dataset/session management.
#' Uses the GLOBAL PGX reactiveValues — loading a dataset in the copilot's
#' left panel propagates to every other board, and selecting a dataset in
#' DataView/Home propagates into the copilot chat.

copilot_parse_features <- function(features) {
  if (is.null(features) || !nzchar(trimws(features))) {
    return(NULL)
  }

  out <- trimws(unlist(strsplit(features, ",")))
  out[nzchar(out)]
}

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

#' FIFO-prune persisted chat sessions to a maximum count.
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

#' Replay persisted user/assistant text turns into a shinychat instance.
#' Tool turns are intentionally skipped per the no-tool-rendering policy.
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

CopilotBoardServer <- function(id, pgx = NULL, pgx_dir = NULL,
                               chat_dir, docs_dir,
                               maxturns = 10,
                               tiers = "copilot-default",
                               is_data_loaded = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    info("[CopilotBoard] server initialised — tiers=",
         paste(tiers, collapse = ","),
         " chat_dir=", chat_dir,
         " docs_dir=", docs_dir)

    ## --- Agent state ---
    copilot <- shiny::reactiveVal(NULL)
    n_turns <- shiny::reactiveVal(0)
    active_dataset_name <- shiny::reactiveVal(NULL)
    busy <- shiny::reactiveVal(FALSE)
    current_tier <- shiny::reactiveVal(tiers[1])
    restoring <- shiny::reactiveVal(FALSE)

    ## --- Persistent chat store (lazy: no dir is created per-session
    ## until the first successful turn triggers session_save) ---
    chat_store <- omicsagentovi::SessionStore(session_dir = chat_dir)

    ## --- Flush the active session on browser close ---
    session$onSessionEnded(function() {
      wrapper <- shiny::isolate(copilot())
      if (!is.null(wrapper) && shiny::isolate(n_turns()) > 0) {
        try(omicsagentovi::session_save(chat_store, wrapper$agent), silent = TRUE)
      }
    })

    ## --- Populate tier selector once on first flush ---
    copilot_tier_label <- function(t) {
      switch(t,
        `copilot-default` = "Balanced (Default)",
        `copilot-fast`    = "Fast",
        `copilot-deep`    = "Deep Think",
        t
      )
    }
    session$onFlushed(function() {
      tier_choices <- setNames(tiers, vapply(tiers, copilot_tier_label, character(1)))
      shiny::updateSelectInput(session, "tier",
                               choices  = tier_choices,
                               selected = tiers[1])
    }, once = TRUE)

    ## --- Reactive wrapper around global pgx for sub-modules ---
    local_pgx <- shiny::reactive({
      if (is.null(pgx$X)) return(NULL)
      snapshot <- shiny::reactiveValuesToList(pgx)
      copilot_as_pgx(snapshot)
    })

    ## --- Sub-module servers ---
    history_refresh <- shiny::reactiveVal(0)
    selected_dataset <- copilot_panel_datasets_server("datasets", pgx_dir = pgx_dir)
    history <- copilot_panel_history_server(
      "history",
      chat_store = chat_store,
      refresh_trigger = history_refresh
    )
    copilot_panel_docs_server("docs", docs_dir = docs_dir)
    evidence <- copilot_panel_evidence_server("evidence", local_pgx = local_pgx)

    ## --- Plot callback for evidence panel ---
    plot_callback <- function(pgx, plot_type, args) {
      evidence$update_plot(copilot_build_plot(pgx, plot_type, args))
    }

    ## --- Greeting (shown once before any dataset) ---
    session$onFlushed(function() {
      mesg <- "Hi, I'm **BigOmics Copilot**! Load a dataset to get started."
      shinychat::chat_append("chat", mesg)
    }, once = TRUE)

    ## --- Create/reset agent whenever the GLOBAL PGX changes ---
    ## This fires both when the user loads data via DataView/Home and when
    ## the copilot's own left panel writes into `pgx` below.
    shiny::observeEvent(list(pgx$name, pgx$X), {
      shiny::req(pgx$name, pgx$X)
      if (isTRUE(restoring())) return()
      dataset_name <- as.character(pgx$name[[1]])
      pgx_snapshot <- copilot_snapshot_pgx(pgx)

      params <- list(
        name  = dataset_name,
        nsamp = omicsai::omicsai_format_num(ncol(pgx$X), digits = 0),
        ngene = omicsai::omicsai_format_num(nrow(pgx$X), digits = 0)
      )

      if (is.null(copilot())) {
        info("[CopilotBoard] first dataset — creating agent for: ", dataset_name)
        agent_wrapper <- copilot_create_agent(
          pgx = pgx_snapshot,
          plot_callback = plot_callback,
          pgx_dir = pgx_dir,
          tier = current_tier()
        )
        model_name <- tryCatch(agent_wrapper$agent@model, error = function(e) "unknown")
        info("[CopilotBoard] agent created — model=", model_name)
        copilot(agent_wrapper)
        active_dataset_name(dataset_name)
        n_turns(0)
        evidence$clear_plot()
        evidence$clear_table()
        shinychat::chat_clear("chat")
        shinychat::chat_append("chat", omicsai::omicsai_substitute_template(
          "Dataset **{{name}}** loaded — {{nsamp}} samples, {{ngene}} genes. What would you like to know?",
          params
        ))
      } else if (!identical(active_dataset_name(), dataset_name)) {
        info("[CopilotBoard] dataset switch in place: ",
             active_dataset_name(), " -> ", dataset_name)
        copilot()$set_pgx(pgx_snapshot)
        active_dataset_name(dataset_name)
        shinychat::chat_append("chat", omicsai::omicsai_substitute_template(
          "_Switched to dataset **{{name}}** — {{nsamp}} samples, {{ngene}} genes._",
          params
        ))
      }
    }, ignoreNULL = TRUE)

    ## --- When user picks a dataset in the copilot's left panel, load it
    ## into the GLOBAL PGX so it propagates to all boards. ---
    shiny::observeEvent(selected_dataset(), {
      dataset_path <- selected_dataset()
      shiny::req(dataset_path)
      shiny::req(file.exists(dataset_path))

      ## Avoid reloading if already the active dataset
      dataset_name <- sub("[.]pgx$", "", basename(dataset_path))
      if (!is.null(pgx$name) && pgx$name == dataset_name) {
        info("[CopilotBoard] dataset already active, skipping reload: ", dataset_name)
        return()
      }

      info("[CopilotBoard] loading PGX into global state: ", dataset_path)
      ## NOTE: this safe-load guard lives here for now because the Copilot
      ## board is the only call-site that loads PGX files outside the normal
      ## Loading wizard. Before merging to devel we will evaluate moving these
      ## checks into playbase so pgx.load / pgx.initialize themselves fail
      ## loudly on corrupt or non-PGX input (currently pgx.initialize silently
      ## returns NULL for wrong-shape lists, and pgx.load surfaces raw R I/O
      ## errors).
      load_result <- shiny::withProgress(message = "Loading dataset for Copilot...", {
        tryCatch({
          p <- playbase::pgx.load(dataset_path)
          p <- playbase::pgx.initialize(p)
          ## Post-condition: pgx.initialize silently returns NULL for
          ## non-PGX lists, so validate shape before trusting it.
          if (is.null(p) || !is.list(p) || is.null(p$X)) {
            stop("file did not yield a valid PGX object (missing expression matrix)")
          }
          p$name <- dataset_name
          list(ok = TRUE, pgx = p)
        }, error = function(e) {
          list(ok = FALSE, err = conditionMessage(e))
        })
      })

      if (!load_result$ok) {
        warning("[CopilotBoard] failed to load PGX: ", dataset_path,
                " :: ", load_result$err)
        shiny::showNotification(
          paste0("Failed to load dataset '", dataset_name, "': ", load_result$err),
          type = "error",
          duration = 8
        )
        return()
      }
      loaded_pgx <- load_result$pgx

      ## Copy into global pgx reactiveValues (same semantics as loading_server.R)
      shiny::isolate(sync_rv_from_list(pgx, loaded_pgx))
      info("[CopilotBoard] PGX loaded into global — samples=", ncol(pgx$X), " genes=", nrow(pgx$X))

      ## Signal the app that a new dataset was loaded so sidebar/tabs refresh
      if (is.function(is_data_loaded)) {
        val <- is_data_loaded()
        is_data_loaded(if (is.null(val)) 1L else val + 1L)
      }

      ## Refresh sidebar menu visibility. The is_data_loaded observer in server.R
      ## also calls this, but the direct call here ensures the JS fires immediately
      ## without waiting for the reactive flush cycle.
      bigdash.showTabsGoToDataView(session)
    })

    ## --- Core ask function ---
    ask_copilot <- function(question, showq = TRUE) {
      agent_wrapper <- copilot()
      if (is.null(agent_wrapper)) {
        info("[CopilotBoard] ask_copilot called but no agent — no dataset loaded")
        shinychat::chat_append("chat", "Please load a dataset first.")
        return(NULL)
      }
      if (n_turns() >= maxturns) {
        info("[CopilotBoard] turn limit reached — maxturns=", maxturns)
        shinychat::chat_append("chat",
          paste("You've reached the", maxturns, "turn limit. Click **+ New chat** in the History panel to start a new conversation."))
        return(NULL)
      }
      if (busy()) {
        shinychat::chat_append("chat", "Please wait for the current response to finish.")
        return(NULL)
      }
      if (showq) {
        msg <- list(role = "user", content = question)
        shinychat::chat_append_message("chat", msg, chunk = FALSE)
      }

      busy(TRUE)
      on.exit(busy(FALSE), add = TRUE)

      failed <- FALSE
      response <- tryCatch(
        agent_wrapper$prompt(question),
        error = function(e) {
          failed <<- TRUE
          info("[CopilotBoard] prompt failed: ", conditionMessage(e))
          paste("Copilot error:", conditionMessage(e))
        }
      )

      shinychat::chat_append("chat", response)
      if (!failed) {
        n_turns(n_turns() + 1)
        tryCatch({
          omicsagentovi::session_save(chat_store, agent_wrapper$agent)
          info("[CopilotBoard] saved session turns=", n_turns())
          tryCatch(
            copilot_prune_sessions(chat_store),
            error = function(e) info("[CopilotBoard] prune failed: ", conditionMessage(e))
          )
          history_refresh(history_refresh() + 1L)
        }, error = function(e) {
          info("[CopilotBoard] session_save failed: ", conditionMessage(e))
        })
      }
    }

    ## --- User chat input ---
    shiny::observeEvent(input$chat_user_input, {
      shiny::req(input$chat_user_input)
      ask_copilot(input$chat_user_input, showq = FALSE)
    })

    ## --- Tier change ---
    ## Changing tier rebuilds the agent on the new (model, reasoning_effort)
    ## and resets the chat. History is not preserved in v1 — keeps the code
    ## simple and matches the existing Reset UX.
    shiny::observeEvent(input$tier, {
      new_tier <- input$tier
      if (is.null(new_tier) || !nzchar(new_tier)) return()
      if (identical(new_tier, current_tier())) return()
      if (busy()) {
        shinychat::chat_append("chat",
          "Please wait for the current response to finish before switching tiers.")
        shiny::updateSelectInput(session, "tier", selected = current_tier())
        return()
      }
      info("[CopilotBoard] tier change — ", current_tier(), " -> ", new_tier)
      current_tier(new_tier)

      if (!is.null(pgx$X)) {
        agent_wrapper <- copilot_create_agent(
          pgx = copilot_snapshot_pgx(pgx),
          plot_callback = plot_callback,
          pgx_dir = pgx_dir,
          tier = new_tier
        )
        copilot(agent_wrapper)
      }
      n_turns(0)
      evidence$clear_plot()
      evidence$clear_table()
      shinychat::chat_clear("chat")
      shinychat::chat_append("chat",
        paste0("Switched to **", copilot_tier_label(new_tier), "** tier. Chat reset."))
      shiny::showNotification(
        paste0("Copilot tier: ", copilot_tier_label(new_tier)),
        type = "message", duration = 3
      )
    }, ignoreInit = TRUE)

    ## --- New chat (from History panel) ---
    shiny::observeEvent(history$new_chat_request(), {
      if (history$new_chat_request() == 0) return()
      if (busy()) {
        shinychat::chat_append("chat", "Please wait for the current response to finish.")
        return()
      }
      ## Flush current session one last time.
      wrapper <- copilot()
      if (!is.null(wrapper) && n_turns() > 0) {
        try(omicsagentovi::session_save(chat_store, wrapper$agent), silent = TRUE)
      }
      info("[CopilotBoard] new chat requested — turns so far=", n_turns())
      if (!is.null(pgx$X)) {
        agent_wrapper <- copilot_create_agent(
          pgx = copilot_snapshot_pgx(pgx),
          plot_callback = plot_callback,
          pgx_dir = pgx_dir,
          tier = current_tier()
        )
        copilot(agent_wrapper)
        active_dataset_name(as.character(pgx$name[[1]]))
      } else {
        copilot(NULL)
        active_dataset_name(NULL)
      }
      n_turns(0)
      evidence$clear_plot()
      evidence$clear_table()
      shinychat::chat_clear("chat")
      shinychat::chat_append("chat", "New conversation started. Ask me anything about your data!")
      history_refresh(history_refresh() + 1L)
    }, ignoreInit = TRUE)

    ## --- Restore session ---
    shiny::observeEvent(history$restore_request(), {
      session_id <- history$restore_request()
      shiny::req(session_id)
      if (busy()) {
        shinychat::chat_append("chat", "Please wait for the current response to finish.")
        return()
      }
      info("[CopilotBoard] restoring session: ", session_id)
      restored <- tryCatch(
        omicsagentovi::ovi_restore(
          session_id = session_id,
          session_dir = chat_store@state$session_dir
        ),
        error = function(e) {
          info("[CopilotBoard] restore failed: ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(restored)) {
        shinychat::chat_append("chat", "Failed to restore session.")
        return()
      }

      ## Wrap restored agent in the same shape copilot_create_agent returns.
      restored_agent <- restored
      wrapper <- list(
        prompt = function(text) omicsagentovi::agent_prompt(restored_agent, text),
        agent = restored_agent,
        tier = current_tier(),
        set_pgx = function(new_pgx) {
          omicsagentovi::ovi_set_pgx(restored_agent, copilot_as_pgx(new_pgx))
        },
        reset = function(...) invisible(NULL)
      )

      restoring(TRUE)
      on.exit(restoring(FALSE), add = TRUE)

      copilot(wrapper)
      evidence$clear_plot()
      evidence$clear_table()
      shinychat::chat_clear("chat")

      turns <- tryCatch(restored_agent@chat$get_turns(), error = function(e) list())
      n_user <- copilot_replay_turns("chat", turns)
      n_turns(n_user)
      active_dataset_name(NULL)  # decoupled from current global pgx
      history$restore_request(NULL)
      info("[CopilotBoard] session restored — user turns replayed=", n_user)
    }, ignoreNULL = TRUE)

    ## --- Delete session ---
    shiny::observeEvent(history$delete_request(), {
      session_id <- history$delete_request()
      shiny::req(session_id)
      info("[CopilotBoard] deleting session: ", session_id)
      try(unlink(
        file.path(chat_store@state$session_dir, session_id),
        recursive = TRUE
      ), silent = TRUE)
      history$delete_request(NULL)
      history_refresh(history_refresh() + 1L)
    }, ignoreNULL = TRUE)

    ## --- Example question buttons ---
    shiny::observeEvent(input$ask_describe, {
      shiny::req(copilot())
      ask_copilot("Describe my experiment and the main comparisons being made")
    })
    shiny::observeEvent(input$ask_findings, {
      shiny::req(copilot())
      ask_copilot("Summarize the main findings of this experiment")
    })
    shiny::observeEvent(input$ask_pathways, {
      shiny::req(copilot())
      ask_copilot("List the top enriched pathways and explain how they relate to the experiment")
    })
    shiny::observeEvent(input$ask_biomarkers, {
      shiny::req(copilot())
      ask_copilot("Show the top candidate biomarkers for this experiment")
    })
    shiny::observeEvent(input$ask_plot, {
      shiny::req(copilot())
      ask_copilot("Show me a PCA plot of the samples")
    })
  })
}
