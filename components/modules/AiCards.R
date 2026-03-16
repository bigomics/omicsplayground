##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# ── LLM Model Profiles ───────────────────────────────────────────────
# Master registry of supported text/diagram models.
# etc/OPTIONS LLM_MODELS selects which of these to offer.
# Only models whose API-key env var is set appear in the dropdown.

LLM_MODEL_PROFILES <- list(
  "google:gemini-2.5-flash" = list(
    label = "Gemini 2.5 Flash", group = "Google",
    env_var = "GEMINI_API_KEY",
    defaults = list(temperature = 0.7, max_tokens = 4096L, timeout_seconds = 90L)
  ),
  "groq:openai/gpt-oss-20b" = list(
    label = "GPT-oss 20B", group = "Groq",
    env_var = "GROQ_API_KEY",
    defaults = list(temperature = 0.7, max_tokens = 4096L, timeout_seconds = 30L)
  ),
  "groq:openai/gpt-oss-120b" = list(
    label = "GPT-oss 120B", group = "Groq",
    env_var = "GROQ_API_KEY",
    defaults = list(temperature = 0.7, max_tokens = 4096L, timeout_seconds = 60L)
  ),
  "xai:grok-4-1-fast-non-reasoning" = list(
    label = "Grok 4 Fast", group = "xAI",
    env_var = "XAI_API_KEY",
    defaults = list(temperature = 0.5, max_tokens = 4096L, timeout_seconds = 60L)
  )
)

# ── Image Model Profiles ─────────────────────────────────────────────
# Image-capable models for infographic generation.
# etc/OPTIONS LLM_IMAGE_MODELS selects which to offer.

IMAGE_MODEL_PROFILES <- list(
  "gemini-3.1-flash-image-preview" = list(
    label = "Gemini 3.1 Flash Image", group = "Google",
    env_var = "GEMINI_API_KEY"
  ),
  "gemini-3.1-pro-image-preview" = list(
    label = "Gemini 3.1 Pro Image", group = "Google",
    env_var = "GEMINI_API_KEY"
  )
)

# ── Mirai daemon for async image generation ────────────────────────
# Start one background worker and pre-load omicsai so the first
# image request doesn't time out waiting for daemon + package startup.
mirai::daemons(1)
mirai::mirai({ library(omicsai); "ready" })  # warm up daemon, prevents coldstart timeout
options(omicsai_image_timeout_s = 90)  # fail fast, retry handles cold starts

#' Get the current AI model from user options
#'
#' Reads the LLM model from session userData (set by the root app server).
#' Uses validate/need so the message propagates through the render chain
#' to the UI when no model is configured.
#'
#' @param parent_session The parent Shiny session passed from the board server
#' @return Model ID string
get_ai_model <- function(parent_session) {
  m <- getUserOption(parent_session, "llm_model")
  shiny::validate(shiny::need(
    !is.null(m) && nzchar(m),
    "No AI model configured. Please enable AI and select a model in Settings."
  ))
  m
}

#' Create omicsai_diagram_config with profile defaults
#' @param model_id Character; key into LLM_MODEL_PROFILES
#' @param system_prompt Character; system prompt for diagram generation
#' @param ... Extra args forwarded to omicsai_diagram_config()
#' @return An omicsai_diagram_config object
make_llm_diagram_config <- function(model_id, system_prompt, ...) {
  d <- LLM_MODEL_PROFILES[[model_id]]$defaults
  omicsai::omicsai_diagram_config(
    model = model_id,
    system_prompt = system_prompt,
    temperature = d$temperature,
    ...
  )
}

#' Return x if non-NULL, else y
#' @param x Primary value
#' @param y Fallback value
#' @return x or y
.aicards_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

#' Translate raw API/HTTP error into a user-friendly string
#' @param msg Character; raw error message from tryCatch
#' @return Character; sanitised message for shiny::validate display
.aicards_friendly_error <- function(msg) {
  if (grepl("503|Service Unavailable|high demand", msg, ignore.case = TRUE)) {
    return(paste0(
      "The AI service is temporarily unavailable due to high demand. ",
      "This is usually short-lived \u2014 please try again in a few minutes."
    ))
  }
  if (grepl("429|Too Many Requests|rate.?limit", msg, ignore.case = TRUE)) {
    return("API rate limit reached. Please wait a moment before generating again.")
  }
  if (grepl("timeout|timed.?out", msg, ignore.case = TRUE)) {
    return(paste0(
      "The request timed out. The service may be under load. ",
      "Please try again."
    ))
  }
  if (grepl("Empty response", msg, ignore.case = TRUE)) {
    return("No response received from the AI model (this can happen during high demand). Please try again.")
  }
  if (grepl("No image data", msg, ignore.case = TRUE)) {
    return("The model returned no image. The prompt may have been blocked or the model is unavailable. Please try again.")
  }
  if (grepl("_API_KEY|API.?key|401|Unauthorized", msg, ignore.case = TRUE)) {
    return("The AI service is not configured or the API key is invalid. Please contact your administrator.")
  }
  # Fallback: strip noisy httr2/ellmer call-stack prefixes and keep the message
  msg <- sub("^Error in [^\n]+:\n\\s*", "", msg)
  paste0("Generation failed: ", msg)
}

#' Coerce a reactive, file path, or string into a reactive returning template text
#' @param template_input Reactive, file path, or template string
#' @return A reactive returning the template string
.aicards_as_reactive_template <- function(template_input) {
  if (shiny::is.reactive(template_input)) {
    return(template_input)
  }

  if (is.character(template_input) && length(template_input) == 1) {
    is_path <- grepl("[/\\\\]", template_input) ||
      grepl("\\.(md|txt|prompt)$", template_input, ignore.case = TRUE)

    if (is_path && file.exists(template_input)) {
      template_content <- omicsai::omicsai_load_template(template_input)
      return(shiny::reactive(template_content))
    }

    return(shiny::reactive(template_input))
  }

  stop("template must be a reactive, file path, or template string")
}

#' Coerce a reactive, config object, or NULL into a reactive returning a config
#' @param config_input Reactive, config object, or NULL (uses default_fn)
#' @param config_class Character; expected S3 class name for validation
#' @param default_fn Function; called to create default config when input is NULL
#' @return A reactive returning the config object
.aicards_as_reactive_config <- function(config_input, config_class, default_fn) {
  if (is.null(config_input)) {
    return(shiny::reactive(default_fn()))
  }
  if (shiny::is.reactive(config_input)) {
    return(config_input)
  }
  if (inherits(config_input, config_class)) {
    return(shiny::reactive(config_input))
  }
  stop("config must be a reactive, ", config_class, " object, or NULL")
}

#' Build Quarto YAML frontmatter for PDF export (lualatex + Lato font)
#' @param title Character; document title
#' @return Character vector of YAML frontmatter lines
.aicards_pdf_frontmatter <- function(title = "AI Report") {
  safe_title <- gsub('"', '\\"', title, fixed = TRUE)
  c(
    "---",
    paste0("title: \"", safe_title, "\""),
    "format:",
    "  pdf:",
    "    pdf-engine: lualatex",
    "    documentclass: article",
    "    papersize: a4",
    "    geometry:",
    "      - left=25mm",
    "      - right=20mm",
    "      - top=25mm",
    "      - bottom=25mm",
    "    mainfont: Lato",
    "---",
    ""
  )
}

#' Render markdown text to PDF via Quarto; returns FALSE on failure
#' @param text Character; markdown content
#' @param file Character; output PDF path
#' @param title Character; document title
#' @return TRUE (invisibly) on success, FALSE on failure
.aicards_markdown_to_pdf <- function(text, file, title = "AI Report") {
  txt <- gsub(intToUtf8(8209), "-", as.character(text), fixed = TRUE)
  txt <- trimws(txt)
  if (!nzchar(txt)) {
    return(FALSE)
  }
  if (!grepl("^\\s*#", txt)) {
    txt <- paste0("# ", title, "\n\n", txt)
  }

  tmpdir <- tempfile(pattern = "aicards_pdf_")
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE), add = TRUE)

  md_file <- file.path(tmpdir, "report.qmd")
  out_file <- file.path(tmpdir, "report.pdf")
  writeLines(
    enc2utf8(c(.aicards_pdf_frontmatter(title), txt)),
    con = md_file,
    useBytes = TRUE
  )

  if (requireNamespace("quarto", quietly = TRUE)) {
    ok <- tryCatch({
      quarto::quarto_render(
        input = md_file,
        output_format = "pdf",
        output_file = basename(out_file),
        quiet = TRUE
      )
      # Some quarto versions write output to the working directory rather than
      # alongside the input file; move it into tmpdir if that happened.
      if (!file.exists(out_file)) {
        alt_out <- file.path(getwd(), basename(out_file))
        if (file.exists(alt_out)) file.rename(alt_out, out_file)
      }
      file.exists(out_file)
    }, error = function(e) FALSE)
    if (!isTRUE(ok)) {
      return(FALSE)
    }
  } else {
    quarto_bin <- Sys.which("quarto")
    if (!nzchar(quarto_bin)) {
      return(FALSE)
    }

    cur_wd <- getwd()
    on.exit(setwd(cur_wd), add = TRUE)
    setwd(tmpdir)

    args <- c("render", basename(md_file), "--to", "pdf", "--output", basename(out_file))
    ret <- tryCatch(
      system2(quarto_bin, args = args, stdout = TRUE, stderr = TRUE),
      error = function(e) structure(conditionMessage(e), status = 1)
    )
    status <- attr(ret, "status")
    if (!is.null(status) && status != 0) {
      return(FALSE)
    }
  }

  if (!file.exists(out_file)) {
    return(FALSE)
  }

  file.copy(out_file, file, overwrite = TRUE)
}

#' Fallback PDF renderer using base R graphics (no Quarto dependency)
#' @param text Character; plain text content
#' @param file Character; output PDF path
#' @param title Character; document title
#' @return NULL (writes PDF as side effect)
.aicards_plain_text_to_pdf <- function(text, file, title = "AI Report") {
  txt <- gsub(intToUtf8(8209), "-", as.character(text), fixed = TRUE)
  lines <- unlist(strsplit(txt, "\n", fixed = TRUE), use.names = FALSE)
  wrapped <- unlist(lapply(lines, strwrap, width = 110), use.names = FALSE)
  if (length(wrapped) == 0) {
    wrapped <- ""
  }

  grDevices::pdf(file = file, width = 8.27, height = 11.69, pointsize = 10)
  on.exit(grDevices::dev.off(), add = TRUE)

  page_size <- 62L # approx lines per A4 page at cex=0.72 with mono font
  chunks <- split(wrapped, ceiling(seq_along(wrapped) / page_size))
  for (i in seq_along(chunks)) {
    graphics::plot.new()
    header <- if (i == 1L) c(title, "") else c(paste0(title, " (continued)"), "")
    page_lines <- c(header, chunks[[i]])
    y <- seq(0.98, 0.02, length.out = length(page_lines))
    graphics::text(0.02, y, labels = page_lines, adj = c(0, 1), family = "mono", cex = 0.72)
  }
}

#' Create a Shiny downloadHandler for markdown-to-PDF export
#' @param text_reactive Reactive returning markdown text
#' @param filename Character; base filename (without .pdf)
#' @param title Character; PDF document title
#' @return A shiny::downloadHandler
.aicards_markdown_pdf_download <- function(text_reactive,
                                           filename = "ai-report",
                                           title = "AI Report") {
  shiny::downloadHandler(
    filename = function() paste0(filename, ".pdf"),
    content = function(file) {
      txt <- text_reactive()
      shiny::req(!is.null(txt), nzchar(trimws(as.character(txt))))

      ok <- .aicards_markdown_to_pdf(txt, file = file, title = title)
      if (!isTRUE(ok)) {
        .aicards_plain_text_to_pdf(txt, file = file, title = title)
      }
    }
  )
}

#' Escape triple backticks so text can safely nest inside a fenced code block
#' @param text Character; raw text (NULL coerced to "")
#' @return Character with ``` replaced by ` ` `
.aicards_safe_fenced_block <- function(text) {
  txt <- .aicards_coalesce(text, "")
  gsub("```", "` ` `", txt, fixed = TRUE)
}

#' Assemble final markdown: AI text + disclaimer + optional prompt debug sections
#' @param text Character; AI-generated text
#' @param show_prompt Logical; if TRUE, append system/user prompt sections
#' @param system_prompt Character or NULL; system prompt to display
#' @param user_prompt Character or NULL; user prompt to display
#' @return Character; complete markdown string
.aicards_text_markdown <- function(text,
                                   show_prompt = FALSE,
                                   system_prompt = NULL,
                                   user_prompt = NULL) {
  out <- paste0(
    text,
    "\n\n*Note: This summary was generated by AI and may contain inaccuracies.*"
  )

  if (isTRUE(show_prompt)) {
    if (!is.null(system_prompt) && nzchar(system_prompt)) {
      out <- paste0(
        out,
        "\n\n---\n\n## System Prompt\n\n```text\n",
        .aicards_safe_fenced_block(system_prompt),
        "\n```"
      )
    }
    if (!is.null(user_prompt) && nzchar(user_prompt)) {
      out <- paste0(
        out,
        "\n\n## User Prompt\n\n```text\n",
        .aicards_safe_fenced_block(user_prompt),
        "\n```"
      )
    }
  }

  out
}

#' AI text card server
#'
#' Synchronous text generation module replacing omicsai_text_server wiring.
#'
#' @param id Shiny module namespace ID
#' @param params_reactive Reactive returning template params list
#' @param template_reactive Reactive or scalar template input
#' @param config_reactive Reactive or scalar omicsai_config
#' @param cache Optional omicsai cache object
#' @param watermark Logical; add watermark in PlotModule
#'
#' @return Reactive returning omicsai_result (or NULL)
AiTextCardServer <- function(id,
                             params_reactive,
                             template_reactive,
                             config_reactive,
                             cache = NULL,
                             watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    module_cache <- .aicards_coalesce(cache, omicsai::omicsai_cache_init("mem"))

    get_template <- .aicards_as_reactive_template(template_reactive)
    get_config <- .aicards_as_reactive_config(
      config_reactive,
      config_class = "omicsai_config",
      default_fn = function() {
        model <- Sys.getenv("OMICS_AI_MODEL", "")
        shiny::validate(shiny::need(
          nzchar(model),
          "No AI model configured. Please select a model in Settings."
        ))
        omicsai::omicsai_config(model = model)
      }
    )

    rv <- shiny::reactiveValues(
      status = "idle",
      result = NULL,
      error = NULL,
      prompt = NULL,
      system_prompt = NULL
    )

    shiny::observeEvent(input$generate, {
      params <- params_reactive()
      shiny::req(params)

      template <- get_template()
      shiny::req(template)

      base_config <- get_config()
      shiny::req(base_config)

      rv$status <- "thinking"
      rv$error <- NULL

      selected_style <- .aicards_coalesce(input$style, "short_summary")
      format_instr <- omicsai::omicsai_instructions(paste0("text/", selected_style))
      config <- base_config
      ## Compose format instructions with board-supplied system prompt
      ## (board's system_prompt contains methods context; do not overwrite it)
      existing_sys <- base_config$system_prompt %||% ""
      config$system_prompt <- if (nzchar(existing_sys)) {
        paste(format_instr, existing_sys, sep = "\n\n")
      } else {
        format_instr
      }

      params$style <- selected_style
      rv$system_prompt <- config$system_prompt
      rv$prompt <- tryCatch(
        omicsai::omicsai_substitute_template(template, params),
        error = function(e) conditionMessage(e)
      )

      tryCatch({
        rv$result <- omicsai::omicsai_gen_text(
          template = template,
          params = params,
          config = config,
          cache = module_cache
        )
        rv$status <- "done"
      }, error = function(e) {
        rv$error <- .aicards_friendly_error(conditionMessage(e))
        rv$status <- "error"
      })
    })

    text_render <- function() {
      if (rv$status == "idle") {
        return(shiny::div(
          class = "text-muted",
          "Click 'Generate!' to create an AI summary."
        ))
      }

      if (rv$status == "thinking") {
        return(shiny::div(
          shiny::icon("spinner", class = "fa-spin"),
          " Generating summary..."
        ))
      }

      if (rv$status == "error") {
        shiny::validate(shiny::need(FALSE, rv$error))
      }

      if (rv$status == "done" && !is.null(rv$result)) {
        txt <- .aicards_text_markdown(
          text = rv$result$text,
          show_prompt = isTRUE(input$show_prompt),
          system_prompt = rv$system_prompt,
          user_prompt = rv$prompt
        )
        txt <- opg_markdown_to_html(txt)
        return(shiny::div(class = "ai-summary-content", shiny::HTML(txt)))
      }

      NULL
    }

    PlotModuleServer(
      "text",
      plotlib = "generic",
      plotlib2 = "generic",
      func = text_render,
      func2 = function() {
        shiny::div(text_render(), style = "font-size:22px;")
      },
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI,
      download.pdf = .aicards_markdown_pdf_download(
        text_reactive = shiny::reactive({
          shiny::req(rv$status == "done", !is.null(rv$result))
          .aicards_text_markdown(
            text = rv$result$text,
            show_prompt = isTRUE(input$show_prompt),
            system_prompt = rv$system_prompt,
            user_prompt = rv$prompt
          )
        }),
        filename = "ai-summary",
        title = "AI Summary"
      ),
      pdf.width = 8,
      pdf.height = 5,
      res = c(75, 100),
      add.watermark = watermark
    )

    result_reactive <- shiny::reactive(rv$result)
    result_reactive
  })
}

#' AI diagram card server
#'
#' Synchronous diagram generation module replacing omicsai_diagram_server wiring.
#'
#' @param id Shiny module namespace ID
#' @param params_reactive Reactive returning template params list
#' @param template_reactive Reactive or scalar template input
#' @param config_reactive Reactive or scalar omicsai_diagram_config
#' @param cache Optional omicsai cache object
#' @param trigger_reactive Optional reactive trigger for external generation events
#' @param style Optional style list passed to omicsai_diagram_render()
#'
#' @return Reactive returning omicsai_diagram_result (or NULL)
AiDiagramCardServer <- function(id,
                                params_reactive,
                                template_reactive,
                                config_reactive,
                                cache = NULL,
                                trigger_reactive = NULL,
                                style = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    module_cache <- .aicards_coalesce(cache, omicsai::omicsai_cache_init("mem"))

    get_template <- .aicards_as_reactive_template(template_reactive)
    get_config <- .aicards_as_reactive_config(
      config_reactive,
      config_class = "omicsai_diagram_config",
      default_fn = function() {
        stop("Diagram config is required and must include system_prompt from build_prompt(... )$system")
      }
    )

    rv <- shiny::reactiveValues(
      result = NULL,
      error = NULL,
      prompt = NULL,
      generate_requested = FALSE
    )

    shiny::observeEvent(input$generate, {
      rv$generate_requested <- TRUE
    })

    if (!is.null(trigger_reactive)) {
      shiny::observeEvent(trigger_reactive(), {
        if (isTRUE(trigger_reactive() > 0)) rv$generate_requested <- TRUE
      }, ignoreInit = TRUE)
    }

    shiny::observe({
      shiny::req(isTRUE(rv$generate_requested))
      rv$error <- NULL

      # params_reactive may be NULL while upstream text is still generating;
      # shiny::req() here is intentional — the observer stays subscribed and
      # re-fires automatically once the text reactive becomes available.
      params <- params_reactive()
      shiny::req(params)

      # template and config are structurally fixed; if they are NULL something
      # is mis-configured — surface an error rather than leaving the spinner up.
      template <- get_template()
      config <- get_config()
      if (is.null(template) || is.null(config)) {
        rv$error <- "Diagram card misconfigured: template or config is NULL."
        rv$generate_requested <- FALSE
        return()
      }

      rv$prompt <- tryCatch({
        user_prompt <- omicsai::omicsai_substitute_template(template, params)
        sys_prompt <- config$system_prompt %||% ""
        if (nzchar(sys_prompt)) {
          paste0(
            "## System Prompt\n\n", sys_prompt,
            "\n\n## User Prompt\n\n", user_prompt
          )
        } else {
          user_prompt
        }
      }, error = function(e) conditionMessage(e))

      tryCatch({
        result <- omicsai::omicsai_gen_diagram(
          template = template,
          params = params,
          config = config,
          cache = module_cache
        )
        rv$result <- omicsai::sanitise_diagram_result(result)
        rv$generate_requested <- FALSE
      }, error = function(e) {
        rv$error <- .aicards_friendly_error(conditionMessage(e))
        rv$generate_requested <- FALSE
      })
    })

    diagram_render <- function() {
      if (!is.null(rv$error)) {
        shiny::validate(shiny::need(FALSE, rv$error))
      }

      if (isTRUE(input$show_prompt) && !is.null(rv$prompt)) {
        shiny::validate(shiny::need(FALSE, paste0("--- Diagram Prompt ---\n\n", rv$prompt)))
      }

      shiny::validate(shiny::need(
        !is.null(rv$result),
        "Diagram generation requires an AI Report. Please generate an AI Report first."
      ))

      omicsai::omicsai_diagram_render(rv$result, style = style, layout = input$layout, spread = input$spread)
    }

    diagram_csv_data <- function() {
      res <- rv$result
      if (is.null(res) || is.null(res$edgelist)) return(NULL)

      nodes_list <- res$edgelist$nodes
      edges_list <- res$edgelist$edges

      if (length(nodes_list) == 0) {
        nodes_df <- data.frame(id = character(), label = character(),
                               type = character(), description = character(),
                               stringsAsFactors = FALSE)
      } else {
        nodes_df <- do.call(rbind, lapply(nodes_list, function(n) {
          data.frame(
            id          = n$id %||% NA_character_,
            label       = n$label %||% NA_character_,
            type        = n$type %||% NA_character_,
            description = n[["function"]] %||% NA_character_,
            stringsAsFactors = FALSE
          )
        }))
      }

      if (length(edges_list) == 0) {
        edges_df <- data.frame(from = character(), to = character(),
                               regulation = character(), label = character(),
                               stringsAsFactors = FALSE)
      } else {
        edges_df <- do.call(rbind, lapply(edges_list, function(e) {
          data.frame(
            from       = e$from %||% NA_character_,
            to         = e$to %||% NA_character_,
            regulation = e$regulation %||% NA_character_,
            label      = e$label %||% NA_character_,
            stringsAsFactors = FALSE
          )
        }))
      }

      list(nodes = nodes_df, edges = edges_df)
    }

    PlotModuleServer(
      "diagram",
      plotlib = "visnetwork",
      func = diagram_render,
      func2 = diagram_render,
      csvFunc = diagram_csv_data,
      pdf.width = 10,
      pdf.height = 8,
      res = c(75, 100)
    )

    result_reactive <- shiny::reactive(rv$result)
    result_reactive
  })
}

#' AI image card server
#'
#' Async image generation module using ExtendedTask + mirai.
#' Runs omicsai_gen_image() in a background worker so the Shiny app
#' stays interactive during the 30-60s image generation.
#'
#' @param id Shiny module namespace ID
#' @param params_reactive Reactive returning template params list
#' @param template_reactive Reactive or scalar template input
#' @param config_reactive Reactive or scalar omicsai_image_config
#' @param cache Optional omicsai cache object
#' @param trigger_reactive Optional reactive trigger for external generation events
#' @param watermark Logical; if TRUE, adds watermark logo to PNG/PDF downloads. Default FALSE.
#'
#' @return Reactive returning omicsai_image_result (or NULL)
AiImageCardServer <- function(id,
                              params_reactive,
                              template_reactive,
                              config_reactive,
                              cache = NULL,
                              trigger_reactive = NULL,
                              watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    module_cache <- .aicards_coalesce(cache, omicsai::omicsai_cache_init("mem"))

    get_template <- .aicards_as_reactive_template(template_reactive)
    get_config <- .aicards_as_reactive_config(
      config_reactive,
      config_class = "omicsai_image_config",
      default_fn = function() {
        stop("Image config is required and must include system_prompt from build_prompt(... )$system")
      }
    )

    rv <- shiny::reactiveValues(
      status = "idle",
      result = NULL,
      error = NULL,
      prompt = NULL,
      generate_requested = FALSE,
      retry = FALSE,
      last_invoke_args = NULL
    )

    # ── ExtendedTask: runs omicsai_gen_image in a background mirai worker ──
    image_task <- shiny::ExtendedTask$new(function(full_prompt, config_list, filename) {
      mirai::mirai(
        {
          tryCatch({
            options(omicsai_image_timeout_s = 90)  # fail fast, retry handles cold starts
            cfg <- omicsai::omicsai_image_config(
              model = config_list$model,
              system_prompt = config_list$system_prompt,
              style = config_list$style,
              n_blocks = config_list$n_blocks,
              aspect_ratio = config_list$aspect_ratio,
              image_size = config_list$image_size
            )
            omicsai::omicsai_gen_image(
              template = full_prompt,
              params = NULL,
              config = cfg,
              filename = filename
            )
          }, error = function(e) e)
        },
        full_prompt = full_prompt,
        config_list = config_list,
        filename = filename
      )
    })

    shiny::observeEvent(input$generate, {
      rv$generate_requested <- TRUE
    })

    if (!is.null(trigger_reactive)) {
      shiny::observeEvent(trigger_reactive(), {
        if (isTRUE(trigger_reactive() > 0)) rv$generate_requested <- TRUE
      }, ignoreInit = TRUE)
    }

    # Validate inputs and launch background task
    shiny::observe({
      shiny::req(isTRUE(rv$generate_requested))
      rv$error <- NULL

      params <- params_reactive()
      shiny::req(params)

      template <- get_template()
      config <- get_config()
      if (is.null(template) || is.null(config)) {
        rv$error <- "Image card misconfigured: template or config is NULL."
        rv$generate_requested <- FALSE
        return()
      }

      # Pre-substitute template in the main thread (params may contain reactives)
      full_prompt <- tryCatch(
        omicsai::omicsai_substitute_template(template, params),
        error = function(e) NULL
      )
      if (is.null(full_prompt)) {
        rv$error <- "Failed to assemble image prompt."
        rv$generate_requested <- FALSE
        return()
      }

      # Store prompt for show_prompt toggle
      sys_prompt <- config$system_prompt %||% ""
      rv$prompt <- if (nzchar(sys_prompt)) {
        paste0("## System Prompt\n\n", sys_prompt,
               "\n\n## User Prompt\n\n", full_prompt)
      } else {
        full_prompt
      }

      # Serialize config as a plain list (S7/S3 objects don't survive mirai transport)
      config_list <- list(
        model = config$model,
        system_prompt = config$system_prompt,
        style = config$style,
        n_blocks = config$n_blocks,
        aspect_ratio = config$aspect_ratio,
        image_size = config$image_size
      )

      filename <- tempfile(fileext = ".png")
      rv$generate_requested <- FALSE
      rv$status <- "running"
      rv$error <- NULL
      rv$result <- NULL
      rv$retry <- FALSE
      rv$last_invoke_args <- list(
        full_prompt = full_prompt,
        config_list = config_list,
        filename = filename
      )
      message(sprintf("[INFO][%s] --- [AI-IMAGE] starting background generation (model: %s)",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), config_list$model))
      image_task$invoke(full_prompt, config_list, filename)
    })

    # Collect result when background task completes
    shiny::observeEvent(image_task$result(), {
      task_result <- image_task$result()
      if (inherits(task_result, "error") || inherits(task_result, "errorValue")) {
        msg <- if (inherits(task_result, "error")) {
          conditionMessage(task_result)
        } else {
          as.character(task_result)
        }
        message(sprintf("[WARN][%s] --- [AI-IMAGE] generation failed: %s",
                        format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
        # HTTP 500: server overloaded — skip retry, show friendly message
        is_500 <- inherits(task_result, "httr2_http_500")
        if (is_500) {
          rv$error <- "Image generation server is temporarily overloaded. Please try again \u2014 click 'Generate!'."
          rv$status <- "error"
          rv$retry <- FALSE
          rv$last_invoke_args <- NULL
          return()
        }
        # Retry once on first failure (timeouts, transient errors)
        if (!isTRUE(rv$retry) && !is.null(rv$last_invoke_args)) {
          message(sprintf("[INFO][%s] --- [AI-IMAGE] retrying (attempt 2 of 2)...",
                          format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
          rv$retry <- TRUE
          args <- rv$last_invoke_args
          image_task$invoke(args$full_prompt, args$config_list, args$filename)
        } else {
          rv$error <- .aicards_friendly_error(msg)
          rv$status <- "error"
          rv$retry <- FALSE
          rv$last_invoke_args <- NULL
        }
      } else {
        message(sprintf("[INFO][%s] --- [AI-IMAGE] generation complete (model: %s)",
                        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                        task_result$metadata$model %||% "unknown"))
        rv$result <- task_result
        rv$status <- "done"
        rv$retry <- FALSE
        rv$last_invoke_args <- NULL
      }
    })

    image_render <- function() {
      if (!is.null(rv$error)) {
        shiny::validate(shiny::need(FALSE, rv$error))
      }

      if (isTRUE(input$show_prompt) && !is.null(rv$prompt)) {
        txt <- paste0("## Image Prompt\n\n```text\n",
                      .aicards_safe_fenced_block(rv$prompt),
                      "\n```")
        return(shiny::div(
          class = "ai-prompt-view",
          style = "overflow-y:auto;max-height:100%;padding:1em;",
          shiny::HTML(opg_markdown_to_html(txt))
        ))
      }

      # Show status while background task is running
      if (rv$status == "running") {
        msg <- if (isTRUE(rv$retry)) {
          "This is taking longer than usual\u2026"
        } else {
          "Generating infographic (this can take 60\u2013120s)\u2026"
        }
        return(shiny::div(
          class = "text-muted",
          style = "display:flex;align-items:center;justify-content:center;height:100%;font-size:1.1em;",
          shiny::icon("spinner", class = "fa-spin me-2"),
          msg
        ))
      }

      shiny::validate(shiny::need(
        !is.null(rv$result),
        "Infographic generation requires an AI Report. Please generate an AI Report first."
      ))

      image_path <- rv$result$path
      shiny::validate(shiny::need(
        !is.null(image_path) && file.exists(image_path),
        "Generated image file not found."
      ))

      ext <- tolower(tools::file_ext(image_path))
      mime <- switch(ext,
        jpg = "image/jpeg",
        jpeg = "image/jpeg",
        gif = "image/gif",
        webp = "image/webp",
        svg = "image/svg+xml",
        "image/png"
      )

      data_uri <- base64enc::dataURI(file = image_path, mime = mime)

      shiny::div(
        class = "ai-image-content",
        style = "width:100%;height:100%;display:flex;align-items:center;justify-content:center;",
        shiny::tags$img(
          src = data_uri,
          alt = "AI-generated infographic",
          style = "max-width:100%;max-height:100%;object-fit:contain;"
        )
      )
    }

    PlotModuleServer(
      "image",
      plotlib = "generic",
      plotlib2 = "generic",
      func = image_render,
      func2 = image_render,
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI,
      pdf.width = 8,
      pdf.height = 5,
      res = c(75, 100),
      add.watermark = watermark
    )

    result_reactive <- shiny::reactive(rv$result)
    result_reactive
  })
}

# ── Microsummary ───────────────────────────────────────────────────
# Lightweight per-tab bullet-point LLM takeaways via omicsai_extract().
# UI layer: MicrosummaryUI() and .microsummary_render_alert() in ui/ui-AiCards.R.

# ellmer type spec for structured extraction — guarantees parsed output.
# Uses type_from_schema() instead of type_object() because ellmer 0.4.0's
# as_json() omits additionalProperties:false, which Groq's API requires.
MICROSUMMARY_TYPE_SPEC <- ellmer::type_from_schema(text = '{
  "type": "object",
  "description": "Microsummary bullet points extracted from analysis data",
  "properties": {
    "bullet1": {
      "type": "string",
      "description": "First key finding — one specific sentence with numbers/names"
    },
    "bullet2": {
      "type": "string",
      "description": "Second key finding — one specific sentence with numbers/names"
    },
    "bullet3": {
      "type": "string",
      "description": "Optional third finding, or empty string if only 2 takeaways"
    }
  },
  "required": ["bullet1", "bullet2", "bullet3"],
  "additionalProperties": false
}')

MICROSUMMARY_USER_TEMPLATE <- omicsai::omicsai_instructions("text/microsummary_data")

#' Microsummary Server
#'
#' Renders a unified bs_alert() box per tab: static description text on top,
#' AI-generated bullet takeaways appended below (or a spinner while loading).
#'
#' @param id Module namespace ID (must match MicrosummaryUI id)
#' @param tab_reactive Reactive returning the current tab name
#' @param params_fn Function(tab) returning a named list with dataset_name
#'   and tab_context, or NULL if data is not ready
#' @param model Character; LLM model identifier for microsummary generation
#' @param cache Optional omicsai cache object; defaults to in-memory cache
#' @param board_name Character; board name for prompt context (e.g. "WGCNA")
#' @param tab_names Character vector of tab names with microsummaries
#' @param static_texts Named list of static HTML strings per tab (the original
#'   bs_alert descriptions). Names must match tab_names.
#' @param invalidate_reactive Optional reactive; when it fires, clears the
#'   in-session results cache
MicrosummaryServer <- function(id,
                               tab_reactive,
                               params_fn,
                               model = "groq:openai/gpt-oss-20b",
                               cache = NULL,
                               board_name = "WGCNA",
                               tab_names = character(0),
                               static_texts = list(),
                               invalidate_reactive = NULL,
                               enabled_reactive = NULL) { # reactive gate: when FALSE, no LLM calls are made and bullets revert to static-only
  shiny::moduleServer(id, function(input, output, session) {
    module_cache <- .aicards_coalesce(cache, omicsai::omicsai_cache_init("mem"))

    # Pre-register all outputs — render static-only alert immediately
    for (tn in tab_names) {
      local({
        tab_name <- tn
        oid <- paste0("micro_", gsub(" ", "_", tab_name))
        output[[oid]] <- shiny::renderUI({
          .microsummary_render_alert(static_texts[[tab_name]], bullets = NULL)
        })
      })
    }

    # Track which tabs have been generated (in-session result cache)
    results <- shiny::reactiveValues()

    # Clear session cache when upstream data changes (e.g. WGCNA recomputation)
    if (!is.null(invalidate_reactive)) {
      shiny::observeEvent(invalidate_reactive(), {
        for (nm in names(results)) results[[nm]] <- NULL
      }, ignoreInit = TRUE)
    }

    # Revert all tabs to static-only when AI is disabled mid-session.
    # NOTE: This intentionally resets the display only, not the results cache.
    # Cached LLM bullets in `results[[tab]]` are preserved so that re-enabling
    # AI restores content instantly without a new LLM call. If a fresh call is
    # ever needed on re-enable, invalidate_reactive should be triggered instead.
    if (!is.null(enabled_reactive)) {
      shiny::observeEvent(enabled_reactive(), {
        if (!isTRUE(enabled_reactive())) {
          for (tn in tab_names) {
            local({
              tab_name <- tn
              oid <- paste0("micro_", gsub(" ", "_", tab_name))
              output[[oid]] <- shiny::renderUI({
                .microsummary_render_alert(static_texts[[tab_name]], bullets = NULL)
              })
            })
          }
        }
      }, ignoreInit = TRUE)
    }

    shiny::observe({
      # Guard: skip LLM call when AI is disabled. Re-enabling lifts this guard
      # and serves cached bullets if available (no forced refresh on re-enable).
      if (!is.null(enabled_reactive) && !isTRUE(enabled_reactive())) return()
      tab <- tab_reactive()
      shiny::req(tab)

      output_id <- paste0("micro_", gsub(" ", "_", tab))
      static_html <- static_texts[[tab]]

      # If already generated for this tab, render cached result and stop
      if (!is.null(shiny::isolate(results[[tab]]))) {
        output[[output_id]] <- shiny::renderUI({
          .microsummary_render_alert(static_html, shiny::isolate(results[[tab]]))
        })
        return()
      }

      # Extract params — params_fn reads reactives, so this observe
      # re-fires when data becomes available
      params <- tryCatch(params_fn(tab), error = function(e) NULL)
      if (is.null(params)) {
        output[[output_id]] <- shiny::renderUI({
          .microsummary_render_alert(static_html, bullets = NULL)
        })
        return()
      }

      # Show loading state — static text + spinner
      output[[output_id]] <- shiny::renderUI({
        .microsummary_render_alert(static_html, bullets = "loading")
      })

      # Build the user prompt
      prompt_params <- list(
        board_name = board_name,
        tab_name = tab,
        dataset_name = params$dataset_name %||% "unknown dataset",
        tab_context = params$tab_context %||% ""
      )
      rendered_prompt <- omicsai::omicsai_substitute_template(
        MICROSUMMARY_USER_TEMPLATE, prompt_params
      )

      config <- omicsai::omicsai_config(
        model = model,
        system_prompt = omicsai::omicsai_instructions("text/microsummary"),
        temperature = 0.3,
        max_tokens = 1024L
      )

      message("[microsummary] generating for tab '", tab, "'...")

      extracted <- tryCatch({
        omicsai::omicsai_extract(rendered_prompt, type_spec = MICROSUMMARY_TYPE_SPEC, config = config, cache = module_cache)$data
      }, error = function(e) {
        message("[microsummary] error for tab '", tab, "': ", conditionMessage(e))
        NULL
      })

      if (is.null(extracted)) {
        output[[output_id]] <- shiny::renderUI({
          .microsummary_render_alert(static_html, bullets = NULL)
        })
        return()
      }

      # Collect non-empty bullets from the structured response
      bullets <- c(extracted$bullet1, extracted$bullet2, extracted$bullet3)
      bullets <- bullets[nzchar(trimws(bullets))]
      bullets <- head(bullets, 3L)
      message("[microsummary] tab '", tab, "' done: ", length(bullets), " bullets")

      results[[tab]] <- bullets

      output[[output_id]] <- shiny::renderUI({
        .microsummary_render_alert(static_html, bullets)
      })
    })

    invisible(NULL)
  })
}
