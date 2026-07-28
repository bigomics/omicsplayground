##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Pre-redesign version saved as launcher_ui_classic() in launcher_ui_classic.R —
## see docs/superpowers/specs/2026-07-17-app-launcher-redesign-design.md.

launcher_ui <- function(id) {
  ns <- shiny::NS(id) ## namespace
  require(bslib)

  apps <- list(
    ## -------------- dashboards -------------
    list(
      input = "launch_playground",
#      icon = "chart-column",
      icon = "tachograph-digital",      
      label = "Omics Playground",
      description = "Play, visualize & discover",
      rgb = c(73,129,178),
      group = "Dashboards"
    ),
    list(
      input = "launch_qsee",
      icon = "circle-check",
      label = "Qsee/Bsee",
      description = "Visual QC & BC analyzer",
      rgb = c(36,176,148),
      group = "Dashboards"
    ),
    list(
      input = "launch_across",
      icon = "shuffle",
      label = "Cross Explorer",
      description = "Query across datasets",
      rgb = c(180,70,70),
      group = "Dashboards"
    ),
    ## ----------------- apps ---------------
    list(
      input = "launch_idconvert",
      icon = "shuffle",
      label = "ID Converter",
      description = "Annotate features",
      rgb = c(140,35,175),
      group = "Apps"
    ),
    list(
      input = "launch_prism",
      icon = "wand-magic-sparkles",
      label = "SmartPrism",
      description = "AI-generated figures",
      rgb = c(190,120,50),
      group = "Apps"
    ),
    list(
      input = "launch_pca",
      icon = "bullseye",
      label = "PCA explorer",
      description = "Interactive PCA",
      rgb = c(190,190,50),
      group = "Apps"
    ),
    list(
      input = "launch_kegg",
      icon = "object-group",
      label = "KEGG Retriever",
      description = "Download KEGG sets",
      rgb = c(190,190,50),
      group = "Apps"
    ),
    list(
      input = "launch_compare",
      icon = "balance-scale",
      label = "PGX Compare",
      description = "Compare PGX side-by-side",
      rgb = c(100, 150, 200),
      group = "Apps"
    )
  )

  app_tile <- function(app) {
    rgb2hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)
    bgcolor <- function(rgb) {
      hue1 <- rgb2hex(rgb[1],rgb[2],rgb[3])
      hue2 <- rgb2hex(0.6*rgb[1],0.6*rgb[2],0.6*rgb[3])
      glue::glue("background: linear-gradient(160deg,{hue1},{hue2});")
    }
    shiny::actionButton(
      ns(app$input),
      label = shiny::tagList(
        shiny::icon(app$icon, class = "app-tile-icon"),
        shiny::div(class = "app-tile-label", app$label),
        shiny::div(class = "app-tile-description", app$description)
      ),
      class = paste0("app-tile ", if (app$group == "Dashboards") "app-tile-wide" else "app-tile-square"),
      style = bgcolor(app$rgb)
    )
  }

  group_section <- function(title, group_apps) {
    grid_class <- paste0("app-launcher-grid app-launcher-grid-", tolower(title))
    shiny::div(
      class = "app-launcher-group",
      style = "margin-bottom: 40px;",
      shiny::div(
        class = "app-launcher-group-header",
        style = "display: flex; align-items: center; gap: 12px; margin: 0 0 16px 0;",
        shiny::h4(title, style = "margin: 0; white-space: nowrap; color: #888;"),
        shiny::div(style = "flex: 1; height: 1px; background: #ddd;")
      ),
      shiny::div(
        class = grid_class,
        lapply(group_apps, app_tile)
      )
    )
  }

  groups <- unique(vapply(apps, function(a) a$group, character(1)))

  ui <- page_fluid(

    ## header with search bar
    bslib::layout_columns(
      style = "text-align: center; padding: 80px 0 50px 0;",
      col_widths = c(-4,4,-4),
      h1("App Launcher"),
      p("Handy apps for your bioinformatics", style="margin-top: -20px;"),
      shiny::textInput(ns("tools_search"), NULL, placeholder="Search...")
    ),

    ## app launcher grid, grouped by section
    shiny::div(
      style = "padding: 35px 15% 0 15%;",
      lapply(groups, function(g) {
        group_section(g, Filter(function(a) a$group == g, apps))
      })
    )
  )

  return(ui)
}
