##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## Fresh mobile-first app launcher UI designed from scratch for Grok.
## Inspired by iOS/Android home screens: large rectangular tappable buttons
## with vibrant colors, prominent icons, clean typography, and responsive grid.
## Starts from the structure of launcher_ui_classic.R but uses completely new
## markup, styling approach, and layout (no reuse of .app-tile, .app-launcher-grid,
## or existing tile helper logic).
##
## Fully responsive: single column on mobile, 2–4 columns on larger screens.
## Large touch targets (≥60px height), rounded corners, subtle shadows, and
## smooth hover/tap feedback.

launcher_ui_grok <- function(id) {
  ns <- shiny::NS(id)
  require(bslib)
  require(shiny)

  # Define the tools with fresh configuration (no reuse of previous list structure)
  tools <- list(
    list(
      id = "runtool_idconvert",
      icon = "arrows-rotate",
      label = "ID Converter",
      desc = "Convert gene/protein IDs & annotate with latest databases",
      color = "#0ea5e9"  # sky blue
    ),
    list(
      id = "runtool_qc",
      icon = "chart-simple",
      label = "Qsee / Bsee",
      desc = "Visual QC, batch effect detection & data diagnostics",
      color = "#10b981"  # emerald green
    ),
    list(
      id = "runtool_prism",
      icon = "wand-magic-sparkles",
      label = "SmartPrism",
      desc = "AI-powered figure generation from natural language",
      color = "#8b5cf6"  # violet
    )
  )

  # Helper to create one large rectangular launcher button
  tool_button <- function(tool) {
    shiny::actionButton(
      inputId = ns(tool$id),
      label = shiny::tagList(
        shiny::tags$div(class = "grok-tool-icon",
          shiny::icon(tool$icon, class = "fa-3x")
        ),
        shiny::tags$div(class = "grok-tool-label", tool$label),
        shiny::tags$div(class = "grok-tool-desc", tool$desc)
      ),
      class = "grok-tool-btn",
      style = paste0("--tool-color: ", tool$color, ";")
    )
  }

  ui <- bslib::page_fluid(
    class = "grok-tools-page",

    ## header with search bar
    bslib::layout_columns(
      style = "text-align: center; padding: 80px 0 50px 0;",
      col_widths = c(-4,4,-4),
      h1("Smart Tools"),
      p("Handy standalone utilities for your bioinformatics", style="margin-top: -20px;"),
      shiny::textInput(ns("tools_search"), NULL, placeholder="Search...") %>%
        shiny::tagAppendAttributes(
          style = "background: #e5e5ea; border-radius: 12px; border: none;
                   padding: 10px 16px; font-size: 16px; width: 100%;",
          class = "ds-search-input"
        )
      
    ),

    # Main launcher area - responsive grid of large rectangular buttons
    shiny::div(
      class = "grok-launcher",
      style = "padding: 2rem 1rem 3rem 1rem; max-width: 1200px; margin: 0 auto;",
      lapply(tools, tool_button)
    ),

    # Inline CSS (self-contained, no dependency on existing .app-tile rules)
    shiny::tags$style(HTML("
      .grok-tools-page {
        background: #f8fafc;
        min-height: 100vh;
      }

      .grok-header h1 {
        color: #0f172a;
      }

      .grok-launcher {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
        gap: 1.25rem;
      }

      /* Large rectangular buttons - iOS/Android launcher style */
      .grok-tool-btn {
        width: 100% !important;
        height: 180px !important;
        border: none !important;
        border-radius: 24px !important;
        background: linear-gradient(145deg, var(--tool-color), color-mix(in srgb, var(--tool-color), #000 15%)) !important;
        color: white !important;
        box-shadow: 0 10px 25px -5px rgb(0 0 0 / 0.15),
                    0 4px 8px -2px rgb(0 0 0 / 0.1) !important;
        transition: all 0.2s cubic-bezier(0.4, 0, 0.2, 1) !important;
        overflow: hidden;
        position: relative;
        padding: 2rem 1.5rem !important;
        text-align: center;
        display: flex !important;
        flex-direction: column;
        align-items: center;
        justify-content: center;
        font-weight: 500;
      }

      .grok-tool-btn:hover,
      .grok-tool-btn:focus {
        transform: translateY(-8px) scale(1.02);
        box-shadow: 0 20px 40px -10px rgb(0 0 0 / 0.25) !important;
        outline: none;
      }

      .grok-tool-btn:active {
        transform: scale(0.98);
      }

      .grok-tool-icon {
        margin-bottom: 1rem;
        opacity: 0.95;
        transition: transform 0.2s ease;
      }

      .grok-tool-btn:hover .grok-tool-icon {
        transform: scale(1.15);
      }

      .grok-tool-label {
        font-size: 1.35rem;
        font-weight: 700;
        margin-bottom: 0.35rem;
        line-height: 1.2;
      }

      .grok-tool-desc {
        font-size: 0.95rem;
        opacity: 0.9;
        line-height: 1.35;
        max-width: 240px;
      }

      /* Mobile optimizations */
      @media (max-width: 640px) {
        .grok-launcher {
          grid-template-columns: 1fr;
          gap: 1rem;
          padding: 0 0.75rem 2rem 0.75rem;
        }

        .grok-tool-btn {
          height: 160px !important;
          border-radius: 20px !important;
          padding: 1.5rem !important;
        }

        .grok-tool-label {
          font-size: 1.25rem;
        }

        .grok-header {
          padding: 1.5rem 1rem 1rem 1rem;
        }

        .grok-header h1 {
          font-size: 2rem;
        }
      }

      /* Extra large screens */
      @media (min-width: 1200px) {
        .grok-launcher {
          grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
        }
      }
    "))
  )

  return(ui)
}
