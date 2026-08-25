##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Launcher page for dashboards, apps and tools.

launcher_ui <- function(id) {
  ns <- shiny::NS(id) ## namespace
  require(bslib)

  apps <- list(
    ## -------------- dashboards -------------
    list(
      input = "launch_playground",
      icon = "chart-column",
      ##icon = "tachograph-digital",      
      label = "Omics Playground",
      description = "Play, see & discover",
      rgb = c(73,129,178),
      group = "Dashboards",
      badge = NULL      
    ),
    list(
      input = "launch_qsee",
      icon = "circle-check",
      label = "Qsee/Bsee",
      description = "Visual QC & BC analyzer",
      rgb = c(36,176,148),
      group = "Dashboards",
      badge = "pro"
    ),
    list(
      input = "launch_across",
      icon = "shuffle",
      ##icon = "layer-group",
      label = "Across Explorer",
      description = "Query across all your datasets",
      rgb = c(180,70,70),
      group = "Dashboards",
      badge = "beta"
    ),
    list(
      input = "launch_methryl",
      icon = "gem",
      label = "Methryl",
      description = "Methylation explorer",
      rgb = c(155,155,165),
      group = "Dashboards",
      badge = "alpha"
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
    )
  )

  ## -------------- quick actions -------------
  quick_actions <- list(
    list(
      input = "load_example",
      label = "Load example",
      icon = "flask",
      class = "btn-outline-primary"
    ),
    list(
      input = "upload_new_data",
      label = "Upload new data",
      icon = "upload",
      class = "btn-outline-primary"
    ),
    list(
      input = "chat_with_obi",
      label = "Chat with Obi",
      icon = "robot",
      class = "btn-outline-primary"
    )
    # Add more quick actions here following the same structure
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
        if (!is.null(app$badge)) shiny::div(class = "app-tile-badge", app$badge),
        shiny::icon(app$icon, class = "app-tile-icon"),
        shiny::div(class = "app-tile-label", app$label),
        shiny::div(class = "app-tile-description", app$description)
      ),
      class = paste0("app-tile ", if (app$group == "Dashboards") "app-tile-wide" else "app-tile-square"),
      style = bgcolor(app$rgb)
    )
  }

  group_section <- function(title, items, render_item, body_class = NULL, body_style = NULL,
                             margin_bottom = "40px") {
    shiny::div(
      class = "app-launcher-group",
      style = paste0("margin-bottom: ", margin_bottom, ";"),
      shiny::div(
        class = "app-launcher-group-header",
        style = "display: flex; align-items: center; gap: 12px; margin: 0 0 16px 0;",
        shiny::h4(title, style = "margin: 0; white-space: nowrap; color: #888;"),
        shiny::div(style = "flex: 1; height: 1px; background: #ddd;")
      ),
      shiny::div(
        class = body_class,
        style = body_style,
        lapply(items, render_item)
      )
    )
  }

  groups <- unique(vapply(apps, function(a) a$group, character(1)))

  quick_action_button <- function(action) {
    shiny::actionButton(
      ns(action$input),
      label = action$label,
      class = paste("launcher-quick-action btn btn-sm",
        action$class %||% "btn-outline-primary"),
      icon = if (!is.null(action$icon)) shiny::icon(action$icon)
      ##style = "font-size: 1rem; padding: 6px 18px;"
    )
  }

  random_greeting <- function() {
    all.hello <- c(
      "Hello", "Salut", "Hola", "Pivet", "Ni hao", "Ciao", "Hi", "Hoi", "Hej",
      "Yassou", "Selam", "Hey", "Hei", "Grutzi", "Bonjour", "Jak się masz",
      "Namaste", "Salam", "Selamat", "Shalom", "Goeiedag", "Yaxshimusiz"
    )
    greetings <- c(
      "What shall we discover today?",
      "What are you working on today?",
      "Where shall we go today?",
      "What would you like to do today?",
      "Ready to dive into your data?",
      "Ready for new discoveries?"
    )
    hh <- paste0(sample(all.hello, 1),"! ")
    gg <- sample(greetings, 1)
    paste0(hh, gg)
  }

  ui <- page_fluid(
    style = "background-color: #f0f9fd; height: 100vh;",

    div(shiny::actionButton(
      ns("logo_click"),
      shiny::tags$img(src = "static/bigomics-logo-small.png", height = "28px"),
      class = "quick-button", style="border: 0px; background-color: transparent;"),
      style="display: flex; justify-content: flex-end; padding: 20px 0 0 0;"),

    ## header with search bar
    bslib::layout_columns(
      style = "text-align: center; margin: -15px 0 20px 0;",
      col_widths = 12,
      div( random_greeting(), id="welcome-text", style="font-size: 32px;")      
      #div("What would you like to discover today?", id = "welcome-subtext")
    ),
    
    ## app launcher grid, grouped by section
    shiny::div(
      style = "padding: 30px 12% 0 12%;",
      
      group_section(
        "Dashboards", Filter(function(a) a$group == "Dashboards", apps), app_tile,
        body_class = "app-launcher-grid app-launcher-grid-dashboards"
      ),

      ## Quick actions section (using same group header style)
      group_section(
        "Quick actions", quick_actions, quick_action_button,
        body_style = "display: flex; gap: 8px; flex-wrap: wrap; padding-left: 4px;",
        margin_bottom = "30px"
      ),

      ##group_section("Apps", Filter(function(a) a$group == "Apps", apps), app_tile,
      ##  body_class = "app-launcher-grid app-launcher-grid-apps"),

    ),

    ## credits
    shiny::div(
      style = "margin-top: 70px; width: 100%;",
      creditsCarousel() 
    )
    
  )

  return(ui)
}
