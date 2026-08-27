##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Tablet-friendly app launcher UI — inspired by mobile home screens but
## scaled for wider screens.  Deliberately different from launcher_ui.R and
## launcher_ui_claude.R.

launcher_ui_ds <- function(id) {
  ns <- shiny::NS(id)
  require(bslib)

  apps <- list(
    list(
      input    = "runtool_idconvert",
      glyph    = "\u21C4",
      label    = "ID Converter",
      desc     = "Convert & annotate features",
      gradient = "linear-gradient(135deg, #007AFF, #5856D6)"
    ),
    list(
      input    = "runtool_qc",
      glyph    = "\u2B50",
      label    = "Qsee/Bsee",
      desc     = "Visual QC & batch effects",
      gradient = "linear-gradient(135deg, #34C759, #00C7BE)"
    ),
    list(
      input    = "runtool_prism",
      glyph    = "\u2728",
      label    = "SmartPrism",
      desc     = "AI-generated figures",
      gradient = "linear-gradient(135deg, #AF52DE, #FF2D55)"
    )
  )

  ## ── App icon (home-screen style, tablet size) ──────────────────────────
  app_icon <- function(app) {
    shiny::actionLink(
      ns(app$input),
      label = NULL,
      style = "text-decoration: none; display: block; cursor: pointer;"
    ) %>%
      shiny::tagAppendChild(
        shiny::div(
          class = "ds-app-icon-wrapper",
          style = "display: flex; flex-direction: column; align-items: center; gap: 10px;",
          shiny::div(
            class = "ds-app-icon-bg",
            style = sprintf(
              "width: 84px; height: 84px; border-radius: 20px; background: %s;
               display: flex; align-items: center; justify-content: center;
               font-size: 36px; color: #fff; box-shadow: 0 6px 20px rgba(0,0,0,0.12);
               transition: transform .15s ease;",
              app$gradient
            ),
            HTML(app$glyph)
          ),
          shiny::div(
            class = "ds-app-label",
            style = "font-size: 14px; font-weight: 500; color: #1c1c1e;
                     text-align: center; line-height: 1.3;",
            app$label
          )
        )
      ) %>%
      shiny::tagAppendAttributes(
        class = "ds-app-cell",
        style = "display: inline-flex; flex-direction: column; align-items: center;"
      )
  }

  ## ── UI ──────────────────────────────────────────────────────────────────
  ui <- bslib::page_fluid(
    style = "background: #f2f2f7; min-height: 100vh; padding: 0 0 60px 0;
             font-family: -apple-system, 'Helvetica Neue', Helvetica, Arial, sans-serif;",

    ## ── Status bar (iPhone-style) ────────────────────────────────────────
    shiny::div(
      class = "ds-status-bar",
      style = "display: flex; justify-content: space-between;
               max-width: 720px; margin: 0 auto;
               padding: 10px 20px 4px 20px; font-size: 12px; font-weight: 600;
               color: #1c1c1e;",
      shiny::div("9:41"),
      shiny::div(
        style = "display: flex; gap: 5px; align-items: center;",
        shiny::icon("signal"),
        shiny::icon("wifi"),
        shiny::icon("battery-full")
      )
    ),

    ## ── Header ───────────────────────────────────────────────────────────
    shiny::div(
      style = "text-align: center; padding: 20px 20px 8px 20px;",
      shiny::div(
        style = "font-size: 32px; font-weight: 700; color: #1c1c1e;",
        "Smart Tools"
      ),
      shiny::div(
        style = "font-size: 15px; color: #8e8e93; margin-top: 4px;",
        "Handy standalone utilities for your bioinformatics"
      )
    ),

    ## ── Search bar ────────────────────────────────────────────────────────
    shiny::div(
      style = "max-width: 480px; margin: 0 auto; padding: 20px 20px 0 20px;",
      shiny::textInput(
        ns("tools_search"),
        label = NULL,
        placeholder = "Search",
        width = "100%"
      ) %>%
        shiny::tagAppendAttributes(
          style = "background: #e5e5ea; border-radius: 12px; border: none;
                   padding: 10px 16px; font-size: 16px; width: 100%;",
          class = "ds-search-input"
        )
    ),

    ## ── App grid ──────────────────────────────────────────────────────────
    shiny::div(
      style = "max-width: 720px; margin: 0 auto; padding: 32px 20px 0 20px;",
      shiny::div(
        style = "display: grid;
                 grid-template-columns: repeat(auto-fill, minmax(120px, 1fr));
                 gap: 28px 16px; justify-items: center;",
        lapply(apps, app_icon)
      )
    ),

    ## ── Hover / active styles ────────────────────────────────────────────
    shiny::tags$head(
      shiny::tags$style(HTML("
        .ds-app-icon-wrapper:hover .ds-app-icon-bg {
          transform: scale(1.08);
        }
        .ds-app-icon-wrapper:active .ds-app-icon-bg {
          transform: scale(0.92);
        }
        .ds-search-input .form-control {
          background: #e5e5ea !important;
          border: none !important;
          border-radius: 12px !important;
          padding: 12px 16px !important;
          font-size: 16px !important;
          box-shadow: none !important;
        }
        .ds-search-input .form-control:focus {
          outline: none !important;
          box-shadow: 0 0 0 2px rgba(0,122,255,0.3) !important;
        }
        .ds-search-input .form-group {
          margin-bottom: 0 !important;
        }
      "))
    )
  )

  return(ui)
}