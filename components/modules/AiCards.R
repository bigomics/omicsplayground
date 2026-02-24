##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

.aicards_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

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

.aicards_markdown_html <- function(text) {
  opg_markdown_to_html(text)
}

.aicards_pdf_frontmatter <- function(title = "AI Report") {
  # pdf-engine: lualatex requires TeX Live with lualatex installed.
  # mainfont: Lato requires the Lato font package (texlive-fonts-extra on Debian/Ubuntu).
  # On systems without these, .aicards_markdown_to_pdf() returns FALSE and the caller
  # falls back to .aicards_plain_text_to_pdf() automatically.
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

.aicards_safe_fenced_block <- function(text) {
  txt <- .aicards_coalesce(text, "")
  gsub("```", "` ` `", txt, fixed = TRUE)
}

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

.aicards_gen_text <- function(template, params, config, cache) {
  omicsai::omicsai_gen_text(
    template = template,
    params = params,
    config = config,
    cache = cache
  )
}

.aicards_gen_diagram <- function(template, params, config, cache) {
  omicsai::omicsai_gen_diagram(
    template = template,
    params = params,
    config = config,
    cache = cache
  )
}

.aicards_gen_infographic <- function(template, params, config, cache) {
  omicsai::omicsai_gen_image(
    template = template,
    params = params,
    config = config,
    cache = cache
  )
}

.aicards_diagram_render <- function(result, style = NULL, layout = NULL, spread = 0.5) {
  omicsai::omicsai_diagram_render(result, style = style, layout = layout, spread = spread)
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
        omicsai::omicsai_config(
          model = Sys.getenv("OMICS_AI_MODEL", "ollama:llama3.2")
        )
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
      rv$status <- "thinking"
      rv$error <- NULL

      params <- params_reactive()
      shiny::req(params)

      template <- get_template()
      shiny::req(template)

      base_config <- get_config()
      shiny::req(base_config)

      selected_style <- .aicards_coalesce(input$style, "short")
      format_instr <- omicsai::omicsai_instructions(paste0("format_", selected_style))
      config <- base_config
      config$system_prompt <- format_instr

      params$style <- selected_style
      rv$system_prompt <- config$system_prompt
      rv$prompt <- tryCatch(
        omicsai::omicsai_substitute_template(template, params),
        error = function(e) conditionMessage(e)
      )

      tryCatch({
        rv$result <- .aicards_gen_text(
          template = template,
          params = params,
          config = config,
          cache = module_cache
        )
        rv$status <- "done"
      }, error = function(e) {
        rv$error <- conditionMessage(e)
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
        return(shiny::div(
          class = "text-danger",
          shiny::icon("exclamation-triangle"),
          " Error: ",
          rv$error
        ))
      }

      if (rv$status == "done" && !is.null(rv$result)) {
        txt <- .aicards_text_markdown(
          text = rv$result$text,
          show_prompt = isTRUE(input$show_prompt),
          system_prompt = rv$system_prompt,
          user_prompt = rv$prompt
        )
        txt <- .aicards_markdown_html(txt)
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
      default_fn = omicsai::omicsai_diagram_config
    )

    rv <- shiny::reactiveValues(
      result = NULL,
      error = NULL,
      generate_requested = FALSE
    )

    request_generation <- function() {
      rv$generate_requested <- TRUE
    }

    shiny::observeEvent(input$generate, {
      request_generation()
    })

    if (!is.null(trigger_reactive)) {
      shiny::observeEvent(trigger_reactive(), {
        request_generation()
      }, ignoreInit = TRUE)
    }

    shiny::observe({
      shiny::req(isTRUE(rv$generate_requested))
      rv$error <- NULL
      params <- params_reactive()
      shiny::req(params)

      template <- get_template()
      shiny::req(template)

      config <- get_config()
      shiny::req(config)

      tryCatch({
        result <- .aicards_gen_diagram(
          template = template,
          params = params,
          config = config,
          cache = module_cache
        )
        rv$result <- omicsai::sanitise_diagram_result(result)
        rv$generate_requested <- FALSE
      }, error = function(e) {
        rv$error <- conditionMessage(e)
        rv$generate_requested <- FALSE
      })
    })

    diagram_render <- function() {
      if (!is.null(rv$error)) {
        shiny::validate(shiny::need(FALSE, paste("Diagram generation failed:", rv$error)))
      }

      shiny::validate(shiny::need(
        !is.null(rv$result),
        "Click 'Generate Diagram' to create an AI diagram."
      ))

      .aicards_diagram_render(rv$result, style = style, layout = input$layout, spread = input$spread)
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
#' Synchronous image generation module replacing omicsai_image_server wiring.
#'
#' @param id Shiny module namespace ID
#' @param params_reactive Reactive returning template params list
#' @param template_reactive Reactive or scalar template input
#' @param config_reactive Reactive or scalar omicsai_image_config
#' @param cache Optional omicsai cache object
#' @param trigger_reactive Optional reactive trigger for external generation events
#'
#' @return Reactive returning omicsai_image_result (or NULL)
AiImageCardServer <- function(id,
                              params_reactive,
                              template_reactive,
                              config_reactive,
                              cache = NULL,
                              trigger_reactive = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    module_cache <- .aicards_coalesce(cache, omicsai::omicsai_cache_init("mem"))

    get_template <- .aicards_as_reactive_template(template_reactive)
    get_config <- .aicards_as_reactive_config(
      config_reactive,
      config_class = "omicsai_image_config",
      default_fn = omicsai::omicsai_image_config
    )

    rv <- shiny::reactiveValues(
      result = NULL,
      error = NULL,
      generate_requested = FALSE
    )

    request_generation <- function() {
      rv$generate_requested <- TRUE
    }

    shiny::observeEvent(input$generate, {
      request_generation()
    })

    if (!is.null(trigger_reactive)) {
      shiny::observeEvent(trigger_reactive(), {
        request_generation()
      }, ignoreInit = TRUE)
    }

    shiny::observe({
      shiny::req(isTRUE(rv$generate_requested))
      rv$error <- NULL
      params <- params_reactive()
      shiny::req(params)

      template <- get_template()
      shiny::req(template)

      config <- get_config()
      shiny::req(config)

      tryCatch({
        rv$result <- .aicards_gen_infographic(
          template = template,
          params = params,
          config = config,
          cache = module_cache
        )
        rv$generate_requested <- FALSE
      }, error = function(e) {
        rv$error <- conditionMessage(e)
        rv$generate_requested <- FALSE
      })
    })

    image_render <- function() {
      if (!is.null(rv$error)) {
        return(shiny::div(
          class = "text-danger",
          shiny::icon("exclamation-triangle"),
          " Error: ",
          rv$error
        ))
      }

      shiny::validate(shiny::need(
        !is.null(rv$result),
        "Click 'Generate Image' to create an AI infographic."
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
      res = c(75, 100)
    )

    result_reactive <- shiny::reactive(rv$result)
    result_reactive
  })
}
