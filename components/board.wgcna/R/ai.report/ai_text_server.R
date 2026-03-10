## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.

# =============================================================================
# AI Text Server for WGCNA
# =============================================================================
# Shiny module that handles text generation for both summary and report modes.
# All pure data extraction lives in ai_report_data.R.
#
# Prompt assembly uses omicsai structured prompt classes (summary_prompt,
# report_prompt) with frag() for deferred template loading.

#' WGCNA AI Text Server
#'
#' Handles text generation for both summary and report modes.
#' In summary mode, generates a single-module summary via omicsai_gen_text().
#' In report mode, generates a full integrated report via omicsai_gen_text()
#' with structured prompt assembly.
#'
#' @param id Module namespace ID
#' @param wgcna Reactive returning WGCNA results object
#' @param pgx PGX object (non-reactive)
#' @param controls List of control reactives (from ai_report_controls_server)
#' @param parent_session Parent Shiny session (for getUserOption)
#'
#' @return List with reactives: text (current mode), report_text (report-only)
wgcna_ai_text_server <- function(id, wgcna, pgx, controls, parent_session,
                                  progress_reactive = NULL) {
  moduleServer(id, function(input, output, session) {

    # ---- Shared: prompt template paths ----

    summary_data_path   <- file.path(BOARD_PROMPTS_DIR, "wgcna_summary_data.md")
    interpretation_path <- file.path(BOARD_PROMPTS_DIR, "wgcna_interpretation.md")

    # ---- Prompt caches (for instant toggle without regeneration) ----
    summary_prompt_cache <- shiny::reactiveVal(NULL)
    report_prompt_cache <- shiny::reactiveVal(NULL)

    # ---- SUMMARY MODE: per-module summary ----

    ai_summary <- shiny::eventReactive(
      list(controls$trigger(), controls$selected_module()),
    {
      mode <- controls$mode()
      if (mode != "summary" || controls$trigger() < 1) return(NULL)

      w <- wgcna()
      module <- controls$selected_module()
      shiny::req(w, module)
      model <- get_ai_model(parent_session)

      style <- controls$summary_style() %||% "short_summary"
      data_params <- wgcna_build_summary_params(w, module, pgx)
      organism <- pgx$organism %||% NULL

      p <- omicsai::summary_prompt(
        role    = omicsai::frag("system_base"),
        task    = omicsai::frag(paste0("text/", style)),
        species = omicsai::omicsai_species_prompt(organism),
        context = omicsai::frag(interpretation_path),
        data    = omicsai::frag(summary_data_path, data_params)
      )
      bp <- omicsai::build_prompt(p)

      ## Cache the prompt for instant toggle
      summary_prompt_cache(paste0(
        "# SYSTEM PROMPT\n\n", bp$system,
        "\n\n---\n\n",
        "# USER MESSAGE\n\n", bp$board
      ))

      cfg <- omicsai::omicsai_config(model = model, system_prompt = bp$system)
      result <- tryCatch(
        omicsai::omicsai_gen_text(bp$board, config = cfg),
        error = function(e) {
          shiny::validate(shiny::need(FALSE, .aicards_friendly_error(conditionMessage(e))))
        }
      )
      result$text
    }, ignoreNULL = FALSE)

    # ---- REPORT MODE: single-shot integrated report ----

    ai_report <- shiny::eventReactive(controls$trigger(), {
      mode <- controls$mode()
      if (mode != "report" || controls$trigger() < 1) return(NULL)

      w <- wgcna()
      shiny::req(w)
      model <- get_ai_model(parent_session)

      # Use coordinator progress if available, otherwise create local
      p <- if (!is.null(progress_reactive)) progress_reactive() else NULL
      local_progress <- is.null(p)
      if (local_progress) {
        p <- shiny::Progress$new()
      }

      ## Step 1: Build structured data tables
      message(sprintf("[INFO][%s] --- [AI-REPORT] extracting module data...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      if (!is.null(p)) p$set(message = "Extracting module data...", value = 0.05)
      tables <- wgcna_build_report_tables(w, pgx)

      ## Step 2: Classify module signal strength
      message(sprintf("[INFO][%s] --- [AI-REPORT] classifying modules...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      if (!is.null(p)) p$set(message = "Classifying modules...", value = 0.10)
      ranking <- wgcna_rank_modules(w)

      ## Step 3: Assemble data content (domain data only, no instructions)
      message(sprintf("[INFO][%s] --- [AI-REPORT] assembling prompt...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      if (!is.null(p)) p$set(message = "Assembling prompt...", value = 0.15)
      data_content <- omicsai::collapse_lines(
        "## Module Signal Classification",
        omicsai::omicsai_format_ranking(ranking),
        "---",
        "## Input Data",
        tables$text,
        sep = "\n\n"
      )

      ## Step 4: Build structured prompt
      board_rules_path <- file.path(BOARD_PROMPTS_DIR, "wgcna_report_rules.md")
      organism <- pgx$organism %||% NULL

      ## Pre-render methods/context (reused in Step 6 as deterministic appendix)
      methods_text <- wgcna_build_methods(w, pgx)

      rp <- omicsai::report_prompt(
        role        = omicsai::frag("system_base"),
        task        = omicsai::frag("text/report"),
        species     = omicsai::omicsai_species_prompt(organism),
        context     = omicsai::frag(interpretation_path),
        board_rules = omicsai::frag(board_rules_path),
        data        = data_content
      )
      bp <- omicsai::build_prompt(rp)

      ## Cache the full prompt for instant toggle
      report_prompt_cache(paste0(
        "# SYSTEM PROMPT\n\n", bp$system,
        "\n\n---\n\n",
        "# USER MESSAGE\n\n", bp$board
      ))

      ## Step 5: Generate report via single LLM call
      ## Use omicsai_gen_text (NOT gen_report — that would double-load text/report.md)
      message(sprintf("[INFO][%s] --- [AI-REPORT] generating AI report (model: %s)...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), model))
      if (!is.null(p)) p$set(message = "Generating AI report...", value = 0.20)

      cfg <- omicsai::omicsai_config(model = model, system_prompt = bp$system, max_tokens = 8192L)
      cache <- omicsai::omicsai_cache_init("mem")
      result <- tryCatch(
        omicsai::omicsai_gen_text(bp$board, config = cfg, cache = cache),
        error = function(e) {
          shiny::validate(shiny::need(FALSE, .aicards_friendly_error(conditionMessage(e))))
        }
      )

      ## Safety net: prepend H1 if missing (same as gen_report does)
      if (!grepl("^# ", result$text)) {
        result$text <- paste0("# Analysis Report\n\n", result$text)
      }

      ## Step 6: Append deterministic methods section
      full_report <- paste(result$text, methods_text, sep = "\n\n")

      message(sprintf("[INFO][%s] --- [AI-REPORT] report text complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      if (!is.null(p)) p$set(message = "Report text complete!", value = 0.65)
      # Don't close progress — coordinator closes after diagram completes
      if (local_progress && !is.null(p)) p$close()
      full_report
    }, ignoreNULL = FALSE)

    # ---- Unified text content (instant toggle between result and prompt) ----

    ai_text <- shiny::reactive({
      mode <- controls$mode()
      show <- isTRUE(controls$show_prompt())

      if (mode == "summary") {
        if (show) summary_prompt_cache() else ai_summary()
      } else {
        if (show) report_prompt_cache() else ai_report()
      }
    })

    list(
      text = ai_text,
      report_text = ai_report
    )
  })
}
