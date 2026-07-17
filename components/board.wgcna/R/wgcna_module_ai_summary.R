##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' WGCNA module AI summary UI
#'
#' Renders an AI summary card for the selected WGCNA module. Used by all four
#' WGCNA tabs (standard, consensus, preservation, multi-omics).
#'
#' @param id Shiny module namespace ID.
#' @param label Optional PlotModule label.
#' @param title Card title.
#' @param info.text Help text.
#' @param caption Caption text.
#' @param height Card height vector c(default, fullscreen).
#' @param width Card width vector c(default, fullscreen).
#'
#' @return Shiny UI tags.
wgcna_module_ai_summary_ui <- function(id,
                                       label = "",
                                       title = "",
                                       info.text = "",
                                       caption = "",
                                       height,
                                       width,
                                       show_save = FALSE) {
  AiTextCardUI(
    id = id,
    title = title,
    caption = caption,
    info.text = info.text,
    height = height,
    width = width,
    show_save = show_save
  )
}

#' WGCNA module AI summary server
#'
#' Serves the per-module AI summary card in two modes:
#'
#' * **Durable** (`variant = "wgcna"` / `"wgcna_mox"`): reads the precomputed
#'   summary from `pgx$ai$<variant>$extras$<module>` and renders it instantly.
#'   "Regenerate" runs a fresh summary and writes it back into the pgx so it
#'   overrides the stored one, mirroring the AI Studio reports.
#' * **Ephemeral** (`variant = NULL`, consensus / preservation): the WGCNA
#'   object is computed live from user input and never persisted, so the summary
#'   is generated on demand and not stored.
#'
#' All prompt assembly is delegated to `playbase::wgcna.module_summary_prompt()`
#' so the board keeps no parallel module-data extraction code.
#'
#' @param id Shiny module namespace ID.
#' @param wgcna Reactive returning the current WGCNA object.
#' @param pgx Current PGX object (reactiveValues).
#' @param r_module Reactive returning the selected module name.
#' @param parent_session Parent Shiny session with AI user options.
#' @param watermark Logical; add watermark in PlotModule.
#' @param user_email Optional user email for telemetry attribution.
#' @param variant Storage slot for durable summaries (`"wgcna"` /
#'   `"wgcna_mox"`), or `NULL` for ephemeral (consensus / preservation).
#' @param board_type WGCNA flavour selecting the prompt/extraction: `"standard"`,
#'   `"multiomics"`, `"consensus"`, or `"preservation"`. Playbase handles each
#'   object's shape (including multi-omics layer resolution) internally.
#'
#' @return Reactive returning the latest omicsai result, or NULL.
wgcna_module_ai_summary_server <- function(id,
                                           wgcna,
                                           pgx,
                                           r_module,
                                           parent_session,
                                           watermark = FALSE,
                                           user_email = NULL,
                                           variant = NULL,
                                           board_type = "standard",
                                           save_pgx = NULL) {
  # The whole WGCNA object is passed to playbase, which selects the right
  # extractor per board_type (multi-omics resolves the owning layer and computes
  # cross-omics coordination from all layers, so it needs the full object).
  wgcna_for_module <- shiny::reactive({
    res <- wgcna()
    module <- r_module()
    shiny::req(res, module)
    list(obj = res, annot = res$annot %||% pgx$genes)
  })

  # Assemble the summary prompt from playbase (single source of truth). The
  # board_type routes consensus / preservation objects to their own extractors.
  prompt_parts <- shiny::reactive({
    wm <- wgcna_for_module()
    playbase::wgcna.module_summary_prompt(
      pgx, wm$obj, r_module(),
      annot = wm$annot, board_type = board_type
    )
  })

  # The card's live/regenerate path uses the already-built board message plus a
  # config carrying the WGCNA system prompt.
  template <- shiny::reactive(prompt_parts()$board)
  config <- shiny::reactive({
    omicsai::omicsai_config(
      model = get_ai_model(parent_session),
      system_prompt = prompt_parts()$system
    )
  })

  # Durable variants: seed the card from the stored summary and persist any
  # regeneration back into the pgx. Ephemeral variants pass NULL for both.
  prefetch_reactive <- NULL
  on_generated <- NULL
  if (!is.null(variant)) {
    prefetch_reactive <- shiny::reactive({
      module <- r_module()
      shiny::req(module)
      pgx$ai[[variant]]$extras[[module]]
    })
    on_generated <- function(result) {
      module <- shiny::isolate(r_module())
      if (is.null(module)) return(invisible(NULL))
      entry <- list(
        summary    = result$text,
        prompt     = shiny::isolate(tryCatch(
          paste0("# SYSTEM\n\n", prompt_parts()$system,
                 "\n\n---\n\n# BOARD\n\n", prompt_parts()$board),
          error = function(e) NULL
        )),
        usage      = result$metadata$usage,
        model      = result$metadata$model %||% NA_character_,
        created_at = as.numeric(Sys.time()),
        edited     = FALSE,
        edited_at  = ""
      )
      if (is.null(pgx$ai)) pgx$ai <- list()
      if (is.null(pgx$ai[[variant]])) pgx$ai[[variant]] <- list()
      if (is.null(pgx$ai[[variant]]$extras)) pgx$ai[[variant]]$extras <- list()
      pgx$ai[[variant]]$extras[[module]] <- entry
    }
  }

  # Persist the in-memory pgx (with any regenerated summaries) to disk, mirroring
  # the AI Studio "Save edits" action. Only durable variants with a save handler
  # get a Save control; ephemeral consensus/preservation cards leave it NULL.
  on_save <- NULL
  if (!is.null(variant) && !is.null(save_pgx)) {
    on_save <- function() save_pgx(pgx)
  }

  # Render the on-demand AI summary card. Telemetry is tagged per variant.
  AiTextCardServer(
    id = id,
    params_reactive = shiny::reactive(list()),
    template_reactive = template,
    config_reactive = config,
    cache = omicsai::omicsai_cache_init("mem"),
    watermark = watermark,
    user_email = user_email,
    telemetry_source = paste0(variant %||% "wgcna", "_summary"),
    # Block generation when AI is off (deployment licence or the Settings
    # "Enable AI" switch, published on session$userData by AppSettingsBoard).
    # Displaying a stored durable summary bypasses this gate (no generation).
    enabled_reactive = shiny::reactive(
      isTRUE(opt$ENABLE_AI) &&
        !isFALSE(getUserOption(parent_session, "ai_enabled"))
    ),
    disabled_message = "Please enable AI to generate module summaries.",
    prefetch_reactive = prefetch_reactive,
    on_generated = on_generated,
    on_save = on_save
  )
}
