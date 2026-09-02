##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

AppSettingsUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  ## IMPORANT: most id's do not use ns() to be accessible at top level
  ## app_server and opg_server.

  ## FORCE_BASIC deploys start basic too, otherwise the full menu renders and
  ## then flips once the login observer runs (server.R). A per-user FORCE_BASIC
  ## is only known after login, so that case still flips.
  initial_is_basic <- (opt$USER_LEVEL == "BASIC" || isTRUE(as.logical(opt$FORCE_BASIC)))

  panel <- bslib::navset_pill_list(
    id = ns("tabs1"),
    widths = c(2,10),
    bslib::nav_panel(
      "App settings",
      bslib::layout_columns(
        height = "100%",
        row_heights = 1,
        ## AI card moved to "AI Features" tab; 1 content card + spacer = c(5,7)
        col_widths = c(4, 8),
        bslib::card(
          bslib::card_header("App settings"),
          bslib::card_body(
            gap = '0.3em',
            bslib::input_switch("enable_beta", "Enable beta features", value=FALSE),
            bslib::input_switch("enable_info", "Show info boxes", value=TRUE),
            bslib::input_switch("menu_basic",  "Basic menu", value=initial_is_basic),
            selector_switch(
              class = "card-footer-checked",
              label = "Show captions",
              is.checked = FALSE
            ),
            withTooltip(
              shiny::selectInput(
                inputId = "selected_labeltype",
                label = "Label type:",
                choices = c("feature", "symbol", "name"),
                selected = "feature",
                width = "100%"
              ),
              "Choose a label type to be displayed in the plots",
              placement = "right", options = list(container = "body")
            )
          )
        ),

        div()
      )
    ),

    ## AI Features tab -------------------------------------------------------
    ## Inputs (wired in appsettings_server.R):
    ##   ns("enable_ai")            — AI on/off toggle (same id as before)
    ##   ns("ai_provider")          — provider dropdown
    ##   ns("ai_api_key")           — API key (collapsed for bigomics)
    ##   ns("ai_base_url")          — base URL (only for custom provider)
    ##   ns("ai_test_load")         — validates BYOK credentials and narrows menus
    ##   ns("ai_test_status")       — status badge for the latest validation
    ##   ns("llm_reports")          — model for AI reports
    ##   ns("llm_images")           — model for image generation
    ##   ns("llm_copilot_deep")     — model for Copilot deep tier
    ##   ns("llm_copilot_balanced") — model for Copilot balanced tier
    bslib::nav_panel(
      "AI Features",
      bslib::layout_columns(
        height = "100%",
        col_widths = c(4, 8),
        bslib::card(
          bslib::card_header("AI Provider"),
          bslib::card_body(
            gap = "0.3em",
            bslib::input_switch(ns("enable_ai"), "Enable AI", value = TRUE),
            shiny::selectInput(
              inputId = ns("ai_provider"),
              label = "AI provider",
              choices = opt$AI_PROVIDERS,
              selected = "bigomics",
              width = "100%"
            ),
            shiny::conditionalPanel(
              condition = "input.ai_provider != 'bigomics'",
              ns = ns,
              shiny::passwordInput(
                inputId = ns("ai_api_key"),
                label = "API key",
                width = "100%"
              )
            ),
            shiny::conditionalPanel(
              condition = "input.ai_provider == 'custom'",
              ns = ns,
              shiny::textInput(
                inputId = ns("ai_base_url"),
                label = "Endpoint base URL",
                width = "100%"
              )
            ),
            shiny::conditionalPanel(
              condition = "input.ai_provider != 'bigomics'",
              ns = ns,
              bslib::layout_columns(
                col_widths = c(7, 5),
                class = "mt-2",
                shiny::actionButton(
                  inputId = ns("ai_test_load"),
                  label = "Test & load models",
                  icon = shiny::icon("plug"),
                  class = "btn-outline-secondary",
                  width = "100%"
                ),
                shiny::uiOutput(ns("ai_test_status"))
              )
            ),

            ## Data-sharing consent. BigOmics-managed backend only: on a BYOK
            ## key the user's own provider account governs their data policy,
            ## and it would be misleading to imply we control it. Off by
            ## default; persisted per user (see components/modules/AiConsent.R)
            ## because a consent that resets every login is not a consent.
            shiny::conditionalPanel(
              condition = "input.ai_provider == 'bigomics'",
              ns = ns,
              shiny::tags$hr(class = "my-2"),
              bslib::input_switch(
                ns("ai_share_data"),
                ## Names the party that actually does the training. BigOmics
                ## trains no models; what the switch changes is which upstream
                ## providers the prompt may be routed to.
                "Allow model providers to train on my Copilot conversations",
                value = FALSE
              ),
              shiny::tags$small(
                class = "text-muted",
                paste(
                  "Off by default: prompts go to providers under a",
                  "no-training, zero-retention policy. Turning this on lets",
                  "them keep and train on your conversations. Your chat",
                  "history is stored in your own workspace either way, and",
                  "usage accounting records token counts only — never",
                  "message content."
                )
              )
            )
          )
        ),

        bslib::card(
          bslib::card_header("AI Models"),
          bslib::card_body(
            gap = "0.3em",
            ## Always show the model menus, even for the BigOmics-managed
            ## backend — but they're disabled (shinyjs, see the lock
            ## observer in appsettings_server.R) for bigomics, since users
            ## don't choose which models BigOmics runs, only see them.
            shiny::conditionalPanel(
              "input.enable_ai",
              ns = ns,
              shiny::selectInput(
                inputId = ns("llm_reports"),
                label = "Reports model",
                choices = opt$AI_MENU_REPORTS,
                selected = 1,
                width = "100%"
              ),
              shiny::selectInput(
                inputId = ns("llm_images"),
                label = "Image AI model",
                choices = opt$AI_MENU_IMAGES,
                selected = 1,
                width = "100%"
              ),
              shiny::selectInput(
                inputId = ns("llm_copilot_deep"),
                label = "Copilot deep model",
                choices = opt$AI_MENU_COPILOT_DEEP,
                selected = 1,
                width = "100%"
              ),
              shiny::selectInput(
                inputId = ns("llm_copilot_balanced"),
                label = "Copilot balanced model",
                choices = opt$AI_MENU_COPILOT_BALANCED,
                selected = 1,
                width = "100%"
              )
            )
          )
        )
      )
    ),

    bslib::nav_panel(
      "New features",
      bslib::layout_columns(
        height = "100%",
        PlotModuleUI(
          ns("newfeatures"),
          outputFunc = htmlOutput,
          info.text = "New features and changes",
          title = "New features"
        )
      )
    ),

    bslib::nav_panel(
      "Package versions",
      bslib::layout_columns(
        height = "100%",
        TableModuleUI(
          ns("packages"),
          info.text = "Packages and versions used in Omics Playground.",
          title = "Package versions"
        )
      )
    ),

    # Plot Colors #####
    ## DEVMODE-only: with the plot editor's body built lazily, theme changes
    ## reach plot modules only through DOM-bound inputs, so a palette change
    ## does not apply until that plot's editor has been opened (PR #1846).
    if (isTRUE(opt$DEVMODE)) bslib::nav_panel(
      "Plot Colors",
      bslib::layout_columns(
        height = "100%",
        row_heights = 1,
        col_widths = c(3,3,3,3),
        bslib::card(
          bslib::card_header("Directional Colors"),
          bslib::card_body(
            gap = "8px",
            colourpicker::colourInput(ns("theme_primary"), "Primary (Up / High)", "#f23451"),
            colourpicker::colourInput(ns("theme_secondary"), "Secondary (Down / Low)", "#3181de"),
            colourpicker::colourInput(ns("theme_neutral"), "Mid / Zero (heatmap)", "#eeeeee")
          )
        ),
        bslib::card(
          bslib::card_header("Chart Elements"),
          bslib::card_body(
            gap = "8px",
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
            gap = "8px",
            colourpicker::colourInput(ns("theme_accent"), "Accent (one significant)", "#e3a45a"),
            colourpicker::colourInput(ns("theme_success"), "Success (both significant)", "#5B9B5B"),
            colourpicker::colourInput(ns("theme_ns_color"), "Not significant", "#eeeeee"),
            colourpicker::colourInput(ns("theme_line"), "Enrichment Line", "#00EE00")
          )
        ),
        bslib::card(
          shiny::actionButton(ns("theme_reset"), "Reset to defaults", icon = shiny::icon("rotate-left"), class = "btn-outline-secondary")
        )
      )
    ),

    bslib::nav_panel(
      "Resource info",
      bslib::layout_columns(
        col_widths = 12,
        height = "100%",
        user_table_resources_ui(ns("resources"))
      )
    )
  )


  OmicsBoardUI(
    ns = ns,
    title = "Settings",
    header_margin = "0px",
    panel
  )
  
}
