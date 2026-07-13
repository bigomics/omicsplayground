##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' AI Studio infographic settings panel
#'
#' Builds the settings sidebar used to select report slots, image style, and
#' trigger infographic regeneration.
#'
#' @param id Shiny module namespace ID.
#' @return Shiny UI tags for the infographic settings panel.
InfographicSettings <- function(id) {
  ns <- shiny::NS(id)
  info <- paste(
    "AI infographics for the current dataset. These images are AI-generated",
    "and can be inaccurate; please always double-check their content."
  )
  shiny::div(
    style = "padding: 10px 15px;",
    h3("Infographics"),
    shiny::br(),
    div(info),
    shiny::br(),
    shiny::checkboxGroupInput(ns("slot"), "Reports:", choices = character(0),
      width = "100%"),
    shiny::selectInput(ns("style"), "Style:", choices = character(0),
      width = "100%"),
    shiny::actionButton(ns("generate"), "Regenerate infographics",
      class = "btn btn-primary", style = "margin-bottom: 6px;")
  )
}

#' AI Studio infographic preview tabs
#'
#' Builds the tabbed preview area for stored AI infographics. Drug report tabs
#' are added dynamically by \code{InfographicServer()}.
#'
#' @param id Shiny module namespace ID.
#' @return A \code{bslib::navset_tab()} UI object.
InfographicUI <- function(id) {
  ns <- shiny::NS(id)
  bslib::navset_tab(
    id = ns("navset"),
    bslib::nav_panel(title = "Summary",
      bslib::card(
        min_height = 600,
        full_screen = TRUE,
        div(shiny::imageOutput(ns("summary"), height = "100%", width = "100%"),
          height = "100%", width = "100%", style = "text-align: center;")
      )
    ),
    bslib::nav_panel(title = "WGCNA",
      bslib::card(
        min_height = 600,
        full_screen = TRUE,
        div(shiny::imageOutput(ns("wgcna"), height = "100%", width = "100%"),
          height = "100%", width = "100%", style = "text-align: center;")
      )
    ),
    bslib::nav_panel(title = "moxWGCNA",
      bslib::card(
        min_height = 600,
        full_screen = TRUE,
        div(shiny::imageOutput(ns("wgcna2"), height = "100%", width = "100%"),
          height = "100%", width = "100%", style = "text-align: center;")
      )
    ),
    bslib::nav_panel(title = "MOFA",
      bslib::card(
        min_height = 600,
        full_screen = TRUE,
        div(shiny::imageOutput(ns("mofa"), height = "100%", width = "100%"),
          height = "100%", width = "100%", style = "text-align: center;")
      )
    ),
    bslib::nav_panel(title = "DE",
      bslib::card(
        min_height = 600,
        full_screen = TRUE,
        div(shiny::imageOutput(ns("de"), height = "100%", width = "100%"),
          height = "100%", width = "100%", style = "text-align: center;")
      )
    ),
    bslib::nav_panel(title = "Enrichment",
      bslib::card(
        min_height = 600,
        full_screen = TRUE,
        div(shiny::imageOutput(ns("enrichment"), height = "100%", width = "100%"),
          height = "100%", width = "100%", style = "text-align: center;")
      )
    )
  )
}

#' Resolve AI report slot names to user-facing infographic labels
#'
#' @param pgx_list PGX object as a list.
#' @param modules Character vector of AI report slot names.
#' @return Character vector of labels matching \code{modules}.
infographic_module_labels <- function(pgx_list, modules) {
  labels <- c(
    combined = "Summary",
    wgcna = "WGCNA",
    wgcna_mox = "moxWGCNA",
    mofa = "MOFA",
    de = "Differential Expression",
    pathways = "Enrichment"
  )
  out <- labels[modules]
  out[is.na(out)] <- vapply(modules[is.na(out)], function(slot) {
    if (startsWith(slot, "drugs_")) {
      ai_report_drug_label(pgx_list, slot)
    } else {
      slot
    }
  }, character(1))
  unname(out)
}

#' Resolve an AI report slot to its Shiny image output ID
#'
#' @param slot AI report slot name.
#' @return Shiny output ID for the infographic preview.
infographic_output_id <- function(slot) {
  switch(slot,
    combined = "summary",
    wgcna = "wgcna",
    wgcna_mox = "wgcna2",
    mofa = "mofa",
    de = "de",
    pathways = "enrichment",
    paste0("preview_", slot)
  )
}

#' Build an omicsai infographic prompt for one AI report
#'
#' @param report AI report markdown/text.
#' @param slot AI report slot name.
#' @param pgx_list PGX object as a list.
#' @param style omicsai image style ID.
#' @return Built omicsai prompt with system and board prompt fields.
infographic_build_image_prompt <- function(report, slot, pgx_list, style) {
  label <- infographic_module_labels(pgx_list, slot)
  organism <- pgx_list$organism
  if (is.null(organism) || !nzchar(organism)) organism <- "human"
  clean <- omicsai::omicsai_strip_report_noise(report)
  species <- omicsai::omicsai_image_species_visual(organism)
  prompt <- omicsai::image_prompt(
    role = omicsai::frag("system_base"),
    task = omicsai::frag("image/infographic",
      params = list(board_name = label)),
    species = species,
    style = omicsai::frag(paste0("image/styles/", style)),
    report = clean
  )
  omicsai::build_prompt(prompt)
}

#' Start one asynchronous infographic image generation job
#'
#' Runs \code{omicsai::omicsai_gen_image()} in a mirai worker and returns a
#' promise resolving to a structured job result.
#'
#' @param board Fully built image prompt body.
#' @param config \code{omicsai_image_config} object.
#' @param slot AI report slot name.
#' @param token Dataset token used to reject stale results.
#' @return A promise resolving to a list with \code{slot}, \code{token},
#'   \code{ok}, and either \code{result} or \code{error}.
infographic_image_promise <- function(board, config, slot, token) {
  promises::as.promise(mirai::mirai(
    {
      options(omicsai_image_timeout_s = 90)
      filename <- tempfile(fileext = ".png")
      tryCatch({
        result <- omicsai::omicsai_gen_image(
          template = board,
          params = NULL,
          config = config,
          filename = filename
        )
        list(slot = slot, token = token, ok = TRUE, result = result)
      }, error = function(e) {
        list(slot = slot, token = token, ok = FALSE,
          error = conditionMessage(e), error_class = class(e))
      })
    },
    board = board,
    config = config,
    slot = slot,
    token = token
  ))
}

#' Normalize an async image job result
#'
#' Coerces malformed promise payloads into the same error-shaped result used by
#' provider failures.
#'
#' @param result Raw promise result.
#' @param slot AI report slot name.
#' @param token Dataset token.
#' @return Structured image job result list.
infographic_job_result <- function(result, slot, token) {
  if (is.list(result) && !is.null(result$slot) &&
      !is.null(result$token) && !is.null(result$ok)) {
    return(result)
  }
  list(
    slot = slot,
    token = token,
    ok = FALSE,
    error = paste(as.character(result), collapse = " ")
  )
}

#' AI Studio infographic server
#'
#' Controls infographic tab visibility, selected report slots, async image
#' generation, retry handling, PGX AI schema storage, and progress messages.
#'
#' @param id Shiny module namespace ID.
#' @param pgx ReactiveValues PGX object.
#' @param save_pgx Optional callback used to persist modified PGX state.
#' @return Shiny module server return value.
InfographicServer <- function(id, pgx, save_pgx = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    tmpdir <- file.path(tempdir(), paste0("ai-infographics-", id))
    dynamic_drug_tabs <- character(0)
    mirai_status <- tryCatch(mirai::status(), error = function(e) NULL)
    if (is.null(mirai_status) || identical(mirai_status$daemons, 0L)) {
      mirai::daemons(2)
    }

    # Snapshot reactive PGX state before passing it into helper functions.
    pgx_snapshot <- function() {
      shiny::reactiveValuesToList(pgx)
    }

    # Render one stored infographic slot into its preview image output.
    render_infographic <- function(slot) {
      renderImage({
        img <- ai_infographic_get(pgx_snapshot(), slot)
        shiny::validate(shiny::need(!is.null(img),
          paste("missing",
            # Use a stable display label for static and dynamic report slots.
            infographic_module_labels(pgx_snapshot(), slot),
            "infographic")))
        ai_infographic_render_value(img, tmpdir, paste0("infographic-", slot))
      }, deleteFile = FALSE)
    }

    # Rebuild dynamic drug-connectivity tabs from the current AI report slots.
    sync_drug_tabs <- function(pgx_list = pgx_snapshot()) {
      for (slot in dynamic_drug_tabs) {
        try(bslib::nav_remove("navset", target = slot), silent = TRUE)
      }
      dynamic_drug_tabs <<- character(0)

      target <- "moxWGCNA"
      for (slot in ai_report_drug_slots(pgx_list)) {
        label <- ai_report_drug_label(pgx_list, slot)
        bslib::nav_insert("navset",
          target = target,
          position = "after",
          nav = bslib::nav_panel(
            title = label,
            value = slot,
            bslib::card(
              min_height = 600,
              full_screen = TRUE,
              div(shiny::imageOutput(ns(infographic_output_id(slot)),
                height = "100%",
                width = "100%"), height = "100%", width = "100%",
                style = "text-align: center;")
            )
          )
        )
        local({
          slot_i <- slot
          # Bind each dynamic tab to the matching stored infographic slot.
          output[[infographic_output_id(slot_i)]] <- render_infographic(slot_i)
        })
        dynamic_drug_tabs <<- c(dynamic_drug_tabs, slot)
        target <- slot
      }
      invisible(dynamic_drug_tabs)
    }

    # Sync visible tabs, checkbox choices, and style choices from PGX state.
    update_nav <- function(pgx_list = pgx_snapshot()) {
      sync_drug_tabs(pgx_list)
      slots <- unique(c(ai_report_slots(pgx_list), ai_infographic_slots(pgx_list)))
      targets <- c(
        combined = "Summary",
        wgcna = "WGCNA",
        wgcna_mox = "moxWGCNA",
        mofa = "MOFA",
        de = "DE",
        pathways = "Enrichment"
      )
      for (slot in names(targets)) {
        if (slot %in% slots) {
          bslib::nav_show("navset", targets[[slot]])
        } else {
          bslib::nav_hide("navset", targets[[slot]])
        }
      }

      choices <- slots[slots %in% ai_report_slots(pgx_list)]
      # Reuse central slot labels so checkbox labels match preview tab labels.
      labels <- infographic_module_labels(pgx_list, choices)
      selected <- shiny::isolate(input$slot)
      selected <- intersect(selected, choices)
      if (!length(selected)) selected <- choices
      shiny::updateCheckboxGroupInput(session, "slot",
        choices = stats::setNames(choices, labels),
        selected = selected)

      styles <- tryCatch(omicsai::omicsai_available_image_styles(),
        error = function(e) "bigomics")
      if (!length(styles)) styles <- "bigomics"
      shiny::updateSelectInput(session, "style",
        choices = styles,
        selected = if ("bigomics" %in% styles) "bigomics" else styles[[1L]])
    }

    image_jobs <- shiny::reactiveValues(
      running = FALSE,
      run_id = 0L,
      token = NULL,
      total = 0L,
      done = 0L,
      failed = 0L,
      progress = NULL,
      changed = FALSE
    )

    # Update the Shiny progress object if a batch is active.
    image_progress <- function(detail = NULL) {
      if (is.null(image_jobs$progress)) return(invisible(NULL))
      value <- if (image_jobs$total > 0L) image_jobs$done / image_jobs$total else 0
      if (is.null(detail)) {
        detail <- paste0(image_jobs$done, " / ", image_jobs$total,
          " infographics complete")
      }
      image_jobs$progress$set(
        message = "Generating AI infographics",
        detail = detail,
        value = value
      )
      invisible(NULL)
    }

    # Close progress, persist changed PGX AI state, and notify the user.
    finish_image_jobs <- function(message = NULL, type = "message") {
      if (!is.null(image_jobs$progress)) {
        image_jobs$progress$close()
        image_jobs$progress <- NULL
      }
      image_jobs$running <- FALSE
      if (isTRUE(image_jobs$changed) && !is.null(save_pgx)) save_pgx(pgx)
      update_nav(pgx_snapshot())
      if (!is.null(message)) {
        shiny::showNotification(message, type = type, session = session)
      }
    }

    # Finish the batch once every selected slot has resolved or failed.
    finish_if_all_image_jobs_done <- function() {
      if (image_jobs$done < image_jobs$total) return(invisible(FALSE))
      if (image_jobs$failed > 0L) {
        finish_image_jobs(
          "Some infographics could not be generated. Please try again.",
          type = "warning"
        )
      } else {
        finish_image_jobs("Infographics ready.")
      }
      invisible(TRUE)
    }

    # Store one job result in pgx$ai and advance batch counters.
    store_image_job_result <- function(slot, result = NULL, error = NULL,
                                       style = NULL) {
      updated_pgx <- ai_infographic_set(
        pgx_snapshot(),
        slot,
        result,
        status = if (is.null(error)) "done" else "error",
        error = error,
        style = style
      )
      pgx$ai <- updated_pgx$ai
      image_jobs$changed <- TRUE
      image_jobs$done <- image_jobs$done + 1L
      if (!is.null(error)) image_jobs$failed <- image_jobs$failed + 1L
      invisible(NULL)
    }

    # Handle fulfilled provider jobs, including provider-level error payloads.
    handle_image_job_result <- function(result, slot, token, img_model,
                                        style, run_id, attempt, cred_fn = NULL) {
      # Normalize malformed promise results before reading fields.
      result <- infographic_job_result(result, slot, token)
      if (!isTRUE(image_jobs$running) ||
          !identical(run_id, image_jobs$run_id)) {
        info("[InfographicServer] image job ignored: slot=", slot)
        return(NULL)
      }
      if (!identical(result$token, ai_report_dataset_token(pgx_snapshot()))) {
        finish_image_jobs(
          "Infographics finished for a previous dataset; results were ignored.",
          type = "warning"
        )
        return(NULL)
      }

      if (isFALSE(result$ok)) {
        msg <- result$error
        info("[InfographicServer] image job failed: slot=", result$slot,
          " attempt=", attempt, " error=", msg)
        if (attempt < 3L) {
          image_progress(paste("Retrying",
            # Use the current PGX labels because dynamic drug names can vary.
            infographic_module_labels(pgx_snapshot(), result$slot),
            "infographic..."))
          launch_image_job(result$slot, img_model, style, run_id,
            attempt = attempt + 1L, cred_fn = cred_fn)
          return(NULL)
        }

        store_image_job_result(result$slot, error = msg, style = style)
        image_progress(paste("Could not generate",
          infographic_module_labels(pgx_snapshot(), result$slot),
          "- continuing with remaining infographics..."))
      } else {
        store_image_job_result(result$slot, result$result, style = style)
        info("[InfographicServer] image job complete: slot=", result$slot,
          " done=", image_jobs$done, "/", image_jobs$total)
        image_progress(paste("Finished",
          infographic_module_labels(pgx_snapshot(), result$slot),
          "-", image_jobs$done, "/", image_jobs$total,
          "infographics complete"))
      }

      finish_if_all_image_jobs_done()
      NULL
    }

    # Handle promise rejections from the mirai/promises layer itself.
    handle_image_job_error <- function(err, slot, label, img_model, style,
                                       run_id, attempt, cred_fn = NULL) {
      if (!isTRUE(image_jobs$running) ||
          !identical(run_id, image_jobs$run_id)) {
        return(NULL)
      }
      info("[InfographicServer] image job failed: slot=", slot,
        " attempt=", attempt, " error=", conditionMessage(err))
      if (attempt < 3L) {
        image_progress(paste("Retrying", label, "infographic..."))
        launch_image_job(slot, img_model, style, run_id,
          attempt = attempt + 1L, cred_fn = cred_fn)
        return(NULL)
      }

      store_image_job_result(slot, error = conditionMessage(err), style = style)
      image_progress(paste("Could not generate", label,
        "- continuing with remaining infographics..."))
      finish_if_all_image_jobs_done()
      NULL
    }

    # Build and launch one report-slot image job.
    launch_image_job <- function(slot, img_model, style, run_id, attempt = 1L,
                                 cred_fn = NULL) {
      pgx_list <- pgx_snapshot()
      entry <- ai_report_get(pgx_list, slot)
      if (is.null(entry)) {
        image_jobs$failed <- image_jobs$failed + 1L
        image_jobs$done <- image_jobs$done + 1L
        return(NULL)
      }
      label <- infographic_module_labels(pgx_list, slot)
      token <- image_jobs$token
      image_progress(paste("Writing", label, "infographic. This can take 60-120s..."))

      # Build the omicsai image prompt in the main Shiny process.
      built <- infographic_build_image_prompt(entry$report, slot, pgx_list, style)
      # credentials is stored inside the config object and serialized into the
      # mirai worker along with the rest of the config — never read from session
      # inside the async worker.
      config <- omicsai::omicsai_image_config(
        model = img_model,
        system_prompt = built$system,
        style = style,
        image_size = "1K",
        credentials = cred_fn
      )

      info("[InfographicServer] image job starting: slot=", slot,
        " model=", img_model, " style=", style, " attempt=", attempt)
      # Run the provider call outside the Shiny event loop.
      promise <- infographic_image_promise(built$board, config, slot, token)

      promises::then(
        promise,
        onFulfilled = function(result) {
          handle_image_job_result(result, slot, token, img_model, style,
            run_id, attempt, cred_fn)
        },
        onRejected = function(err) {
          handle_image_job_error(err, slot, label, img_model, style,
            run_id, attempt, cred_fn)
        }
      )
      invisible(promise)
    }

    # Start an async batch for the selected report slots.
    start_image_jobs <- function(slots, img_model, style, cred_fn = NULL) {
      slots <- unique(as.character(slots))
      slots <- slots[!is.na(slots) & nzchar(slots)]
      if (!length(slots)) {
        shiny::showNotification("No AI reports available for infographic generation.",
          type = "warning", session = session)
        return(NULL)
      }

      image_jobs$run_id <- image_jobs$run_id + 1L
      run_id <- image_jobs$run_id
      image_jobs$running <- TRUE
      image_jobs$token <- ai_report_dataset_token(pgx_snapshot())
      image_jobs$total <- length(slots)
      image_jobs$done <- 0L
      image_jobs$failed <- 0L
      image_jobs$changed <- FALSE
      image_jobs$progress <- shiny::Progress$new(session, min = 0, max = 1)
      image_progress("Writing your infographic. This can take 60-120s...")
      for (slot in slots) {
        launch_image_job(slot, img_model, style, run_id, cred_fn = cred_fn)
      }
    }

    output$summary <- render_infographic("combined")
    output$wgcna <- render_infographic("wgcna")
    output$wgcna2 <- render_infographic("wgcna_mox")
    output$mofa <- render_infographic("mofa")
    output$de <- render_infographic("de")
    output$enrichment <- render_infographic("pathways")

    shiny::observeEvent(list(pgx$X, pgx$name, pgx$ai), {
      shiny::req(pgx$X, pgx$name)
      update_nav(pgx_snapshot())
    }, ignoreNULL = FALSE)

    shiny::observeEvent(input$generate, {
      shiny::req(pgx$X, pgx$name)
      if (isTRUE(image_jobs$running)) {
        shiny::showNotification("Infographic generation is already running.",
          type = "message", session = session)
        return(NULL)
      }

      img_model <- getUserOption(session, "img_model")
      if (is.null(img_model) || !nzchar(img_model)) {
        shiny::showNotification("Please enable an image model in user settings.",
          type = "warning", session = session)
        return(NULL)
      }
      cred_fn <- get_ai_credentials(session)

      pgx_list <- pgx_snapshot()
      slots <- input$slot
      if (is.null(slots)) slots <- character(0)
      slots <- as.character(slots)
      slots <- slots[slots %in% ai_report_slots(pgx_list)]
      style <- input$style
      if (is.null(style) || !nzchar(style)) style <- "bigomics"
      start_image_jobs(slots, img_model, style, cred_fn)
    }, ignoreInit = TRUE)

    shiny::observeEvent(pgx$name, {
      if (isTRUE(image_jobs$running)) {
        finish_image_jobs(
          "Infographic generation was cancelled because the dataset changed.",
          type = "warning"
        )
      }
    }, ignoreInit = TRUE)
  })
}
