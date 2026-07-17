##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

AI_TEXT_BUDGETS <- list(
  short_summary = 2048L,
  long_summary = 4096L
)

#' Get the current AI model from user options
#'
#' @param parent_session Parent Shiny session.
#' @return Provider-qualified model ID.
get_ai_model <- function(parent_session) {
  model <- getUserOption(parent_session, "llm_model")
  shiny::validate(shiny::need(
    !is.null(model) && nzchar(model),
    "No AI model configured. Please enable AI and select a model in Settings."
  ))
  model
}

#' Get the current AI credentials from user options
#'
#' @param parent_session Parent Shiny session.
#' @return A nullary credential closure \code{function() key}, or \code{NULL}
#'   when the provider is BigOmics (env-var path) or no key has been set.
get_ai_credentials <- function(parent_session) {
  getUserOption(parent_session, "ai_credentials")
}

#' AI text card server
#'
#' Runs on-demand text generation from a prompt template and reactive params.
#'
#' @param id Shiny module namespace ID.
#' @param params_reactive Reactive returning template params.
#' @param template_reactive Reactive returning a prompt template.
#' @param config_reactive Reactive returning an `omicsai_config`.
#' @param cache Optional omicsai cache object.
#' @param watermark Logical; add watermark in PlotModule.
#' @param user_email Optional user email for telemetry attribution.
#' @param telemetry_source Telemetry `source` tag for this card's generations.
#'
#' @return Reactive returning the latest omicsai result, or NULL.
AiTextCardServer <- function(id,
                             params_reactive,
                             template_reactive,
                             config_reactive,
                             cache = NULL,
                             watermark = FALSE,
                             user_email = NULL,
                             telemetry_source = "card",
                             enabled_reactive = NULL,
                             disabled_message = "AI features are disabled.",
                             prefetch_reactive = NULL,
                             on_generated = NULL,
                             on_save = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    # Prepare card-local state and cache for one on-demand text result.
    module_cache <- if (is.null(cache)) omicsai::omicsai_cache_init("mem") else cache
    rv <- shiny::reactiveValues(
      status = "idle",
      result = NULL,
      error = NULL,
      prompt = NULL,
      system_prompt = NULL
    )

    # Durable-first display: when a precomputed summary exists for the current
    # selection, render it instantly (no LLM call, no AI-enabled gate). Falling
    # back to "idle" when none exists lets the user Generate on demand. Skipped
    # while a live generation is running so a click is never overwritten.
    if (!is.null(prefetch_reactive)) {
      shiny::observeEvent(prefetch_reactive(), ignoreNULL = FALSE, {
        if (identical(rv$status, "running")) return(NULL)
        entry <- prefetch_reactive()
        if (!is.null(entry) && !is.null(entry$summary)) {
          rv$status <- "done"
          rv$result <- list(text = entry$summary)
          rv$error <- NULL
        } else {
          rv$status <- "idle"
          rv$result <- NULL
          rv$error <- NULL
        }
      })
    }

    # Persist the current summaries to disk (durable variants only). The board
    # supplies on_save = function() save_pgx(pgx); we just report success/failure.
    if (!is.null(on_save)) {
      shiny::observeEvent(input$save, {
        ok <- tryCatch(isTRUE(on_save()), error = function(e) FALSE)
        shiny::showNotification(
          if (ok) "Summary saved to dataset." else "Could not save summary.",
          type = if (ok) "message" else "error",
          session = session
        )
      })
    }

    shiny::observeEvent(input$generate, {
      # Guard: block generation when AI is disabled (deployment licence or the
      # user's "Enable AI" switch). Show a friendly message and stop early.
      if (!is.null(enabled_reactive) && !isTRUE(enabled_reactive())) {
        rv$status <- "idle"
        shiny::showNotification(disabled_message, type = "message",
                                session = session)
        return(NULL)
      }
      # Read reactive inputs at click time so tab changes do not auto-run.
      params <- params_reactive()
      template <- template_reactive()
      config <- config_reactive()
      style <- input$style

      # Merge user style instructions into the board-provided system prompt.
      text_config <- .aicards_text_config(config, style)
      style <- text_config$style
      config <- text_config$config

      # Reset result state and store prompt text for optional prompt display.
      rv$status <- "running"
      rv$result <- NULL
      rv$error <- NULL
      rv$prompt <- tryCatch(
        omicsai::omicsai_substitute_template(template, params),
        error = function(e) template
      )
      rv$system_prompt <- tryCatch(
        omicsai::omicsai_config_get(config, "system_prompt"),
        error = function(e) NULL
      )
      model <- tryCatch(omicsai::omicsai_config_get(config, "model"), error = function(e) "unknown")
      started <- Sys.time()
      info("[AiTextCardServer] text generation starting: id=", id,
           " model=", model, " style=", style)

      # Generate in-session; this is intentionally serial for fast text models.
      result <- tryCatch(
        omicsai::omicsai_gen_text(
          template = template,
          params = params,
          config = config,
          cache = module_cache
        ),
        error = function(e) e
      )

      # Convert provider errors into card state instead of throwing in Shiny.
      if (inherits(result, "error")) {
        rv$status <- "error"
        rv$error <- .aicards_friendly_error(conditionMessage(result))
      } else {
        rv$status <- "done"
        rv$result <- result

        # Telemetry must never break the card: usage is omicsai's object/list,
        # passed straight through (no field remapping).
        tryCatch(
          ai_telemetry_record(
            source     = telemetry_source,
            session_id = session$token,
            user_email = user_email %||% NA_character_,
            model      = model,
            usage      = result$metadata$usage
          ),
          error = function(e) NULL
        )

        # Durable variants persist the regenerated summary back into the pgx so
        # it overrides the stored one (AI-Studio style); never break the card.
        if (!is.null(on_generated)) {
          tryCatch(on_generated(result), error = function(e) NULL)
        }
      }
      info("[AiTextCardServer] text generation finished: id=", id,
           " status=", rv$status,
           " elapsed_s=", round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1))
    })

    # Render status, result text, and optional prompt using one body builder.
    contents_text <- shiny::reactive({
      .aicards_text_content(
        status = rv$status,
        result = rv$result,
        error = rv$error,
        prompt = rv$prompt,
        show_prompt = isTRUE(input$show_prompt)
      )
    })

    # Keep PlotModule naming convention used across the application.
    text_RENDER <- function() {
      shiny::div(class = "gene-info", shiny::HTML(opg_markdown_to_html(contents_text())))
    }

    text_RENDER2 <- function() {
      shiny::div(
        shiny::HTML(opg_markdown_to_html(contents_text())),
        style = "font-size:22px;"
      )
    }

    # Delegate card rendering and fullscreen rendering to the existing PlotModule.
    PlotModuleServer(
      "text",
      plotlib = "generic",
      plotlib2 = "generic",
      func = text_RENDER,
      func2 = text_RENDER2,
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI,
      download.fmt = NULL,
      pdf.width = 8,
      pdf.height = 5,
      res = c(75, 100),
      add.watermark = watermark
    )

    shiny::reactive(rv$result)
  })
}

# Compose the style-specific omicsai config for one text generation request.
# This helps us keep prompt policy and token budget logic out of the Shiny event.
.aicards_text_config <- function(config, style) {
  if (is.null(style) || !style %in% names(AI_TEXT_BUDGETS)) {
    style <- "short_summary"
  }

  style_prompt <- tryCatch(
    omicsai::omicsai_instructions(paste0("text/", style)),
    error = function(e) {
      if (identical(style, "long_summary")) {
        return("Write a structured 500-600 word summary with clear biological interpretation.")
      }
      "Write a concise 250-350 word summary focused on the strongest biological signals."
    }
  )
  existing_sys <- tryCatch(
    omicsai::omicsai_config_get(config, "system_prompt"),
    error = function(e) ""
  )
  system_prompt <- if (nzchar(existing_sys)) {
    paste(style_prompt, existing_sys, sep = "\n\n")
  } else {
    style_prompt
  }

  list(
    style = style,
    config = omicsai::omicsai_config(
      base = config,
      system_prompt = system_prompt,
      max_tokens = AI_TEXT_BUDGETS[[style]]
    )
  )
}

# Build the markdown body shown inside the text card.
# Separating the code here makes status rendering independent from API calls.
.aicards_text_content <- function(status,
                                  result,
                                  error,
                                  prompt,
                                  show_prompt = FALSE) {
  if (identical(status, "idle")) {
    return("Select a module and click **Generate** to create an AI summary.")
  }
  if (identical(status, "running")) {
    return("Generating summary. This can take 30-90 seconds...")
  }
  if (identical(status, "error")) {
    return(paste0("**Summary generation failed.**\n\n", error))
  }

  text <- paste0(
    result$text,
    "\n\n*Note: This summary was generated by AI and may contain inaccuracies.*"
  )
  if (isTRUE(show_prompt)) {
    text <- paste0(
      text,
      "\n\n---\n\n## Prompt\n\n```text\n",
      gsub("```", "` ` `", prompt, fixed = TRUE),
      "\n```"
    )
  }
  text
}

# Translate raw provider errors into messages suitable for the card.
# This helps us keep provider-specific error strings away from the UI surface.
.aicards_friendly_error <- function(msg) {
  if (grepl("503|Service Unavailable|high demand", msg, ignore.case = TRUE)) {
    return("The AI service is temporarily unavailable due to high demand. Please try again in a few minutes.")
  }
  if (grepl("429|Too Many Requests|rate.?limit", msg, ignore.case = TRUE)) {
    return("API rate limit reached. Please wait a moment before generating again.")
  }
  if (grepl("timeout|timed.?out", msg, ignore.case = TRUE)) {
    return("The request timed out. Please try again.")
  }
  if (grepl("_API_KEY|API.?key|401|Unauthorized", msg, ignore.case = TRUE)) {
    return("The AI service is not configured or the API key is invalid. Please contact your administrator.")
  }
  paste0("Generation failed: ", sub("^Error in [^\n]+:\\s*", "", msg))
}
