##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# ── Model Choice Builders ─────────────────────────────────────────────
# Build grouped selectInput choices from profile registries
# (LLM_MODEL_PROFILES / IMAGE_MODEL_PROFILES defined in modules/AiCards.R).
# Accepts either a semicolon-separated string (from etc/OPTIONS) or a
# character vector. Only models whose API-key env var is set are offered.

.build_model_choices <- function(profiles, selected_ids = NULL) {
  if (is.character(selected_ids) && length(selected_ids) == 1 && grepl(";", selected_ids)) {
    selected_ids <- trimws(strsplit(selected_ids, ";")[[1]])
    selected_ids <- selected_ids[nzchar(selected_ids)]
  }
  if (!is.null(selected_ids) && length(selected_ids)) {
    profiles <- profiles[intersect(selected_ids, names(profiles))]
  }
  available <- Filter(function(p) nzchar(Sys.getenv(p$env_var)), profiles)
  if (!length(available)) return(character(0))
  split(
    setNames(names(available), vapply(available, `[[`, "", "label")),
    vapply(available, `[[`, "", "group")
  )
}

llm_model_choices <- function(selected_ids = NULL) {
  .build_model_choices(LLM_MODEL_PROFILES, selected_ids)
}

image_model_choices <- function(selected_ids = NULL) {
  .build_model_choices(IMAGE_MODEL_PROFILES, selected_ids)
}

#' AI text card UI
#'
#' PlotModule wrapper for text-only AI output.
#'
#' @param id Shiny module namespace ID
#' @param title Card title
#' @param caption Caption text
#' @param info.text Help text
#' @param height Card height vector c(default, fullscreen)
#' @param width Card width vector c(default, fullscreen)
#'
#' @return Shiny tagList
AiTextCardUI <- function(id,
                         title = "AI Summary",
                         caption = "AI-generated summary.",
                         info.text = "",
                         height = c("100%", TABLE_HEIGHT_MODAL),
                         width = c("auto", "100%")) {
  ns <- shiny::NS(id)

  opts <- shiny::tagList(
    shiny::radioButtons(
      ns("style"),
      "AI summary:",
      choices = c("Short" = "short_summary", "Long" = "long_summary"),
      selected = "short_summary",
      inline = TRUE
    ),
    shiny::checkboxInput(
      ns("show_prompt"),
      "Show prompt",
      FALSE
    ),
    shiny::actionButton(
      ns("generate"),
      "Generate!",
      icon = shiny::icon("refresh"),
      class = "btn-outline-primary"
    )
  )

  PlotModuleUI(
    ns("text"),
    outputFunc = shiny::htmlOutput,
    title = title,
    label = "",
    info.text = info.text,
    options = opts,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("pdf")
  )
}

#' AI diagram card UI
#'
#' PlotModule wrapper for visNetwork AI diagrams.
#'
#' @param id Shiny module namespace ID
#' @param title Card title
#' @param caption Caption text
#' @param info.text Help text
#' @param height Card height vector c(default, fullscreen)
#' @param width Card width vector c(default, fullscreen)
#'
#' @return Shiny tagList
AiDiagramCardUI <- function(id,
                            title = "AI Diagram",
                            caption = "AI-generated diagram.",
                            info.text = "",
                            height = c("100%", TABLE_HEIGHT_MODAL),
                            width = c("auto", "100%")) {
  ns <- shiny::NS(id)

  opts <- shiny::tagList(
    shiny::checkboxInput(
      ns("show_prompt"),
      "Show prompt",
      FALSE
    ),
    shiny::actionButton(
      ns("generate"),
      "Generate Diagram",
      icon = shiny::icon("refresh"),
      class = "btn-outline-primary"
    ),
    shiny::selectInput(
      ns("layout"),
      label = "Layout",
      choices = c(
        "Hierarchical" = "hierarchical",
        "Force-directed (FR)" = "fr",
        "Spring (KK)" = "kk",
        "Auto" = "nicely"
      ),
      selected = "nicely",
      width = "100%"
    ),
    shiny::conditionalPanel(
      condition = sprintf("input['%s'] !== 'hierarchical'", ns("layout")),
      shiny::sliderInput(
        ns("spread"),
        label = "Node spread",
        min = 0, max = 1, value = 0.5, step = 0.05,
        width = "100%"
      )
    )
  )

  PlotModuleUI(
    ns("diagram"),
    plotlib = "visnetwork",
    title = title,
    label = "",
    info.text = info.text,
    options = opts,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "html", "csv")
  )
}

#' AI image card UI
#'
#' PlotModule wrapper for infographic image output.
#'
#' @param id Shiny module namespace ID
#' @param title Card title
#' @param caption Caption text
#' @param info.text Help text
#' @param height Card height vector c(default, fullscreen)
#' @param width Card width vector c(default, fullscreen)
#'
#' @return Shiny tagList
AiImageCardUI <- function(id,
                          title = "AI Image",
                          caption = "AI-generated infographic.",
                          info.text = "",
                          height = c("100%", TABLE_HEIGHT_MODAL),
                          width = c("auto", "100%")) {
  ns <- shiny::NS(id)

  opts <- shiny::tagList(
    shiny::checkboxInput(
      ns("show_prompt"),
      "Show prompt",
      FALSE
    ),
    shiny::actionButton(
      ns("generate"),
      "Generate Image",
      icon = shiny::icon("refresh"),
      class = "btn-outline-primary"
    )
  )

  PlotModuleUI(
    ns("image"),
    outputFunc = shiny::htmlOutput,
    title = title,
    label = "",
    info.text = info.text,
    options = opts,
    caption = caption,
    height = height,
    width = width
  )
}

# ── Microsummary UI ─────────────────────────────────────────────────
# Presentation layer for AI microsummary module.
# MicrosummaryUI()             — Shiny module UI (one uiOutput per tab)
# .microsummary_render_alert() — stateless HTML builder (static text + bullets)

#' Microsummary UI
#'
#' Creates one uiOutput per tab. Each output renders the complete alert box
#' (static description + AI bullets) as a single visual unit.
#'
#' @param id Module namespace ID
#' @param tab_names Character vector of tab names to create outputs for
#' @return A named list of uiOutput tags (one per tab)
MicrosummaryUI <- function(id, tab_names) {
  ns <- shiny::NS(id)
  outputs <- lapply(tab_names, function(tab) {
    shiny::uiOutput(ns(paste0("micro_", gsub(" ", "_", tab))))
  })
  names(outputs) <- tab_names
  outputs
}

#' Render unified alert box: static description + optional AI bullets
#'
#' @param static_html HTML string for the static tab description (or NULL)
#' @param bullets Character vector of bullet strings, "loading" for spinner,
#'   or NULL for static-only
#' @return A bs_alert-style Shiny tag
#' @keywords internal
.microsummary_render_alert <- function(static_html, bullets = NULL) {
  # Build the static part
  parts <- list()
  if (!is.null(static_html) && nzchar(static_html)) {
    parts <- list(shiny::HTML(static_html))
  }

  # Append bullets or spinner
  if (identical(bullets, "loading")) {
    parts <- c(parts, list(
      shiny::tags$hr(style = "margin: 8px 0; border-color: rgba(0,0,0,0.1);"),
      shiny::div(
        style = "font-size: 0.85em; color: #224164;",
        shiny::icon("spinner", class = "fa-spin fa-xs"),
        " Summarizing..."
      )
    ))
  } else if (!is.null(bullets) && length(bullets) > 0) {
    li_html <- paste0("<li>", htmltools::htmlEscape(bullets), "</li>", collapse = "\n")
    parts <- c(parts, list(
      shiny::tags$hr(style = "margin: 8px 0; border-color: rgba(0,0,0,0.1);"),
      shiny::HTML(paste0(
        "<ul style='margin:0; padding-left:1.2em; font-size:0.85em;'>",
        li_html,
        "</ul>"
      ))
    ))
  }

  # Render as a single bs_alert
  bs_alert(
    shiny::tagList(parts),
    conditional = TRUE,
    style = "primary",
    closable = TRUE,
    translate = FALSE,
    html = FALSE
  )
}
