##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Pre-redesign version saved as launcher_ui_classic() in launcher_ui_classic.R —
## see docs/superpowers/specs/2026-07-17-app-launcher-redesign-design.md.

launcher_ui_claude <- function(id) {
  ns <- shiny::NS(id) ## namespace
  require(bslib)

  apps <- list(
    list(
      input = "runtool_playground",
      icon = "gauge-high",
      label = "Omics Playground",
      description = "Play, visualize & discover",
      color = "blue"
    ),
    list(
      input = "runtool_idconvert",
      icon = "shuffle",
      label = "ID Converter",
      description = "Convert & annotate features",
      color = "purple"
    ),
    list(
      input = "runtool_qc",
      icon = "chart-column",
      label = "Qsee/Bsee",
      description = "Visual QC & batch effects",
      color = "teal"
    ),
    list(
      input = "runtool_prism",
      icon = "wand-magic-sparkles",
      label = "SmartPrism",
      description = "AI-generated figures",
      color = "orange"
    )
  )

  app_tile <- function(app) {
    shiny::actionButton(
      ns(app$input),
      label = shiny::tagList(
        shiny::icon(app$icon, class = "app-tile-icon"),
        shiny::div(class = "app-tile-label", app$label),
        shiny::div(class = "app-tile-description", app$description)
      ),
      class = paste0("app-tile app-tile-", app$color)
    )
  }

  ui <- page_fluid(

    ## header with search bar
    bslib::layout_columns(
      style = "text-align: center; padding: 80px 0 50px 0;",
      col_widths = c(-4,4,-4),
      h1("Smart Tools"),
      p("Handy standalone utilities for your bioinformatics", style="margin-top: -20px;"),
      shiny::textInput(ns("tools_search"), NULL, placeholder="Search...")

      ## %>% shiny::tagAppendAttributes(
      ##   style = "background: #e5e5ea; border-radius: 12px; border: none;
      ##              padding: 10px 16px; font-size: 16px; width: 100%;",
      ##   class = "ds-search-input"
      ## )
      
    ),

    ## app launcher grid
    shiny::div(
      class = "app-launcher-grid",
      style = "padding: 0 15%;",
      lapply(apps, app_tile)
    )
  )

  return(ui)
}
