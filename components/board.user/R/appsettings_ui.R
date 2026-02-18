##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

AppSettingsInputs <- function(id) {
  ns <- shiny::NS(id)
  bigdash::tabSettings()
}

AppSettingsUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  div(
    boardHeader(title = "App Settings", info_link = ns("board_info")),
    shiny::tabsetPanel(
      id = ns("tabs1"),
      shiny::tabPanel(
        "App Settings & News",
        bslib::layout_columns(
          height = "calc(100vh - 183px)",
          col_widths = c(6, 6),
          PlotModuleUI(
            ns("newfeatures"),
            outputFunc = htmlOutput,
            info.text = "New features and changes",
            title = "New features"
          ),
          TableModuleUI(
            ns("packages"),
            info.text = "Packages and versions used in Omics Playground.",
            title = "Package versions"
          )
        )
      ),
      # Plot Colors #####
      shiny::tabPanel(
        "Plot Colors",
        bslib::layout_columns(
          height = "calc(100vh - 183px)",
          col_widths = c(4, 4, 4),
          bslib::card(
            bslib::card_header("Directional Colors"),
            bslib::card_body(
              colourpicker::colourInput(ns("theme_primary"), "Primary (Up / High)", "#f23451"),
              colourpicker::colourInput(ns("theme_secondary"), "Secondary (Down / Low)", "#3181de")
            )
          ),
          bslib::card(
            bslib::card_header("Chart Elements"),
            bslib::card_body(
              colourpicker::colourInput(ns("theme_neutral"), "Neutral / Mid", "#eeeeee"),
              colourpicker::colourInput(ns("theme_bar_color"), "Bar Color", "#A6CEE3"),
              shiny::selectInput(
                ns("theme_palette"), "Default Palette",
                choices = c(
                  "muted_light", "default", "light", "dark",
                  "super_light", "super_dark", "muted", "expanded",
                  "highlight_blue", "highlight_red", "highlight_orange",
                  "custom_gradient"
                ),
                selected = "default"
              ),
              shiny::conditionalPanel(
                condition = paste0("input['", ns("theme_palette"), "'] == 'custom_gradient'"),
                colourpicker::colourInput(ns("theme_palette_c1"), "Gradient start", "#3181de"),
                colourpicker::colourInput(ns("theme_palette_c2"), "Gradient middle", "#eeeeee"),
                colourpicker::colourInput(ns("theme_palette_c3"), "Gradient end", "#f23451")
              )
            )
          ),
          bslib::card(
            bslib::card_header("Accent Colors"),
            bslib::card_body(
              colourpicker::colourInput(ns("theme_accent"), "Accent (one significant)", "#e3a45a"),
              colourpicker::colourInput(ns("theme_success"), "Success (both significant)", "#5B9B5B"),
              colourpicker::colourInput(ns("theme_line"), "Enrichment Line", "#00EE00"),
              shiny::hr(),
              shiny::actionButton(ns("theme_reset"), "Reset to defaults", icon = shiny::icon("rotate-left"), class = "btn-outline-secondary")
            )
          )
        )
      ),
      # Resource info #####
      shiny::tabPanel(
        "Resource info",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 183px)",
          user_table_resources_ui(ns("resources"))
        )
      )
    )
  )
}
