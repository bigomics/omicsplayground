##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

AppSettingsBoard <- function(id, auth, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    dbg("[AppSettingsBoard] >>> initializing User Settings...")

    ## module for system resources
    user_table_resources_server("resources", pgx = pgx)

    shiny::observeEvent(input$board_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>User Profile</strong>"),
        shiny::HTML(
          "The User Settings page allows you to change overall settings
                that will alter how the app looks and functions."
        ),
        easyClose = TRUE, size = "l"
      ))
    })

    ## ----------------------------------------------------------------
    ## AI Features — provider, credentials and model selection
    ## ----------------------------------------------------------------

    ## Menu input id -> (catalog key, union-allowlist fallback). The union
    ## fallback is only used when a provider is missing from the catalog.
    ai_menu_choices <- list(
      llm_reports          = list(key = "reports",          union = opt$AI_MENU_REPORTS),
      llm_images           = list(key = "images",           union = opt$AI_MENU_IMAGES),
      llm_copilot_deep     = list(key = "copilot_deep",     union = opt$AI_MENU_COPILOT_DEEP),
      llm_copilot_balanced = list(key = "copilot_balanced", union = opt$AI_MENU_COPILOT_BALANCED)
    )

    provider_menu_choices <- function(provider, input_id, live_models = NULL) {
      spec <- ai_menu_choices[[input_id]]
      if (identical(provider, "custom") && !is.null(live_models) &&
          !identical(spec$key, "images")) {
        return(unique(live_models[!is.na(live_models) & nzchar(live_models)]))
      }

      catalog <- if (is.null(provider) || !nzchar(provider)) {
        NULL
      } else {
        opt$AI_MODELS[[provider]]
      }
      choices <- if (is.null(catalog)) spec$union else catalog[[spec$key]]
      choices <- choices %||% character(0)

      if (!is.null(live_models)) {
        live_choices <- intersect(choices, live_models)
        if (length(live_choices)) {
          return(live_choices)
        }
      }
      choices
    }

    update_ai_model_menus <- function(provider, live_models = NULL) {
      for (input_id in names(ai_menu_choices)) {
        choices <- provider_menu_choices(provider, input_id, live_models)
        shiny::updateSelectInput(
          session, input_id,
          choices = choices,
          selected = if (length(choices)) choices[[1]] else character(0)
        )
      }
    }

    ai_test_status <- shiny::reactiveVal(list(state = "idle", label = "Not tested"))

    output$ai_test_status <- shiny::renderUI({
      status <- ai_test_status()
      state <- status$state %||% "idle"
      label <- status$label %||% "Not tested"
      class <- switch(state,
        ok = "badge bg-success",
        error = "badge bg-danger",
        "badge bg-secondary"
      )
      icon <- switch(state,
        ok = "check-circle",
        error = "exclamation-triangle",
        "circle"
      )
      shiny::tags$div(
        class = "d-flex align-items-center h-100",
        shiny::tags$span(
          class = class,
          style = "width: 100%; padding: 0.55rem 0.6rem;",
          shiny::icon(icon),
          shiny::tags$span(style = "margin-left: 0.35rem;", label)
        )
      )
    })

    ## Warn the user when they switch AI on.
    shiny::observeEvent(input$enable_ai, {
      model <- input$llm_reports
      if (isTRUE(input$enable_ai)) {
        if (is.null(model) || model == "") {
          shinyalert::shinyalert(
            "ERROR",
            "No LLM server available. Please check your settings."
          )
          return(NULL)
        }
        shinyalert::shinyalert("WARNING",
          "Using LLM might expose some of your data to external LLM servers.",
          closeOnClickOutside = TRUE
          # showCancelButton = TRUE
        )
      }
    })
    
    shiny::observeEvent(
      {
        list(
          input$enable_ai,
          input$llm_reports,
          input$llm_images,
          input$llm_copilot_deep,
          input$llm_copilot_balanced
        )
      },
      {
        ## Legacy user-option names ("llm_model"/"img_model") are still read by
        ## the loading + upload boards' report/image generation, so keep writing
        ## them alongside the copilot-tier keys.
        ## Re-check the deployment licence (opt$ENABLE_AI), not just the switch:
        ## on an unlicensed deployment the switch is only greyed (shinyjs::disable),
        ## so its value stays TRUE and would otherwise seed a real llm_model that
        ## the loading/upload auto-report paths then act on.
        if (isTRUE(opt$ENABLE_AI) && isTRUE(input$enable_ai)) {
          setUserOption(session, "llm_model", input$llm_reports)
          setUserOption(session, "img_model", input$llm_images)
          setUserOption(session, "llm_copilot_deep", input$llm_copilot_deep)
          setUserOption(session, "llm_copilot_balanced", input$llm_copilot_balanced)
        } else {
          setUserOption(session, "llm_model", "")
          setUserOption(session, "img_model", "")
          setUserOption(session, "llm_copilot_deep", "")
          setUserOption(session, "llm_copilot_balanced", "")
        }
      }
    )

    ## Store the selected provider and build the session-only credential
    ## closure. It captures the key *value* (not the reactive) and lives only
    ## in session$userData — never logged, never persisted.
    shiny::observeEvent(
      {
        list(input$ai_provider, input$ai_api_key, input$ai_base_url)
      },
      {
        provider <- input$ai_provider
        if (is.null(provider) || !nzchar(provider)) {
          return(NULL)
        }
        setUserOption(session, "ai_provider", provider)

        ## BigOmics uses the env-var path downstream (no key). Any other
        ## provider stores a session-only closure capturing the key *value*
        ## — never the reactive; NULL when the key box is empty.
        key <- if (identical(provider, "bigomics")) NULL else input$ai_api_key
        cred <- if (is.null(key) || !nzchar(key)) {
          NULL
        } else {
          local({ k <- key; function() k })
        }
        setUserOption(session, "ai_credentials", cred)

        if (identical(provider, "custom")) {
          setUserOption(session, "ai_base_url", input$ai_base_url)
        }
      }
    )

    ## Repopulate the four model menus with the selected provider's catalog.
    shiny::observeEvent(input$ai_provider, {
      provider <- input$ai_provider
      if (is.null(provider) || !nzchar(provider)) {
        return(NULL)
      }
      ai_test_status(list(state = "idle", label = "Not tested"))
      update_ai_model_menus(provider)
    })

    ## Test credentials by asking the provider for live model ids. OPG still
    ## owns the menu allowlist/order, so live ids only narrow the JSON-policy
    ## choices; empty/error responses keep the static catalog-derived menus.
    shiny::observeEvent(input$ai_test_load, {
      provider <- input$ai_provider
      if (is.null(provider) || !nzchar(provider) ||
          identical(provider, "bigomics")) {
        ai_test_status(list(state = "idle", label = "Not tested"))
        return(NULL)
      }

      live_models <- omicsai::ai.list_provider_models(
        provider = provider,
        key = input$ai_api_key %||% "",
        base_url = input$ai_base_url %||% NULL
      )
      live_models <- unique(live_models[!is.na(live_models) & nzchar(live_models)])
      matched_models <- unique(unlist(lapply(names(ai_menu_choices), function(input_id) {
        intersect(provider_menu_choices(provider, input_id, live_models), live_models)
      }), use.names = FALSE))
      matched_models <- matched_models[!is.na(matched_models) & nzchar(matched_models)]

      ## Both failure modes leave the menus on the static catalog and flag the
      ## status badge; only the message differs.
      ai_test_fail <- function(title, text) {
        ai_test_status(list(state = "error", label = "Error"))
        update_ai_model_menus(provider)
        shinyalert::shinyalert(title = title, text = text, type = "error",
                               closeOnClickOutside = TRUE)
        NULL
      }

      if (!length(live_models)) {
        return(ai_test_fail(
          "Could not load provider models",
          paste(
            "No models were returned for the selected provider.",
            "Check the API key, provider, and endpoint base URL, then try again.",
            "The model menus were left on the default catalog choices."
          )
        ))
      }

      if (!length(matched_models)) {
        return(ai_test_fail(
          "Provider models do not match Omics Playground",
          paste(
            "The provider returned models, but none match the enabled",
            "Omics Playground menus for this provider.",
            "The model menus were left on the default catalog choices."
          )
        ))
      }

      update_ai_model_menus(provider, live_models = live_models)
      ai_test_status(list(state = "ok", label = "OK"))
      shinyalert::shinyalert(
        title = "API key is correctly set",
        text = sprintf(
          "Loaded %d live model%s from %s. %d model%s match Omics Playground menus.",
          length(live_models),
          if (length(live_models) == 1L) "" else "s",
          provider,
          length(matched_models),
          if (length(matched_models) == 1L) "" else "s"
        ),
        type = "success",
        closeOnClickOutside = TRUE
      )
    })

    ## Admin lock: grey the provider dropdown and the "Enable AI" switch for
    ## non-admins when the deployment pins AI config (AI_PROVIDER_LOCKED), and
    ## always grey them when the deployment is not licensed for AI
    ## (opt$ENABLE_AI == FALSE) — so admins keep control over AI usage.
    shiny::observe({
      locked <- (isTRUE(opt$AI_PROVIDER_LOCKED) && !isTRUE(auth$ADMIN)) ||
        !isTRUE(opt$ENABLE_AI)
      for (input_id in c("ai_provider", "enable_ai")) {
        if (locked) shinyjs::disable(input_id) else shinyjs::enable(input_id)
      }
    })

    ## Effective AI-enabled state = deployment licence (opt$ENABLE_AI) AND the
    ## user/admin "Enable AI" switch. Publish it on session$userData (read by
    ## AI feature guards, e.g. the WGCNA module summary). The AI Studio +
    ## Copilot tabs live in the root "app-sidebar" navset, so the app-level
    ## server show/hides them via the `enable_ai` reactive returned below
    ## (shiny::hideTab namespaces the id, so it must run in the root session).
    shiny::observe({
      ai_on <- isTRUE(opt$ENABLE_AI) && isTRUE(input$enable_ai)
      setUserOption(session, "ai_enabled", ai_on)
    })


    ## ----------------------------------------------------------------
    ## Plot Colors theme handling
    ## ----------------------------------------------------------------
    theme <- get_color_theme()

    ## Load persisted theme when user logs in
    shiny::observeEvent(auth$logged, {
      if (!isTRUE(auth$logged)) {
        return()
      }
      saved <- load_color_theme(auth$user_dir)
      if (is.null(saved)) {
        return()
      }
      ## Apply known keys only (guards against stale/corrupt files)
      for (key in intersect(names(saved), names(COLOR_THEME_DEFAULTS))) {
        theme[[key]] <- saved[[key]]
      }
      ## Sync UI pickers
      colourpicker::updateColourInput(session, "theme_primary", value = theme$primary)
      colourpicker::updateColourInput(session, "theme_secondary", value = theme$secondary)
      colourpicker::updateColourInput(session, "theme_neutral", value = theme$neutral)
      colourpicker::updateColourInput(session, "theme_ns_color", value = theme$ns_color)
      colourpicker::updateColourInput(session, "theme_bar_color", value = theme$bar_color)
      colourpicker::updateColourInput(session, "theme_accent", value = theme$accent)
      colourpicker::updateColourInput(session, "theme_success", value = theme$success)
      colourpicker::updateColourInput(session, "theme_line", value = theme$line)
      shiny::updateSelectInput(session, "theme_palette", selected = theme$palette)
      colourpicker::updateColourInput(session, "theme_palette_c1", value = theme$palette_c1)
      colourpicker::updateColourInput(session, "theme_palette_c2", value = theme$palette_c2)
      colourpicker::updateColourInput(session, "theme_palette_c3", value = theme$palette_c3)
    })

    ## Clear the session credential on logout so a user's key never bleeds
    ## into a later session that reuses this Shiny process.
    shiny::observeEvent(auth$logged, {
      if (isTRUE(auth$logged)) {
        return()
      }
      setUserOption(session, "ai_credentials", NULL)
      setUserOption(session, "ai_provider", "bigomics")
      shiny::updateSelectInput(session, "ai_provider", selected = "bigomics")
    }, ignoreInit = TRUE)

    ## Map UI input IDs to theme keys
    theme_inputs <- list(
      theme_primary   = "primary",
      theme_secondary = "secondary",
      theme_neutral   = "neutral",
      theme_ns_color  = "ns_color",
      theme_bar_color = "bar_color",
      theme_accent    = "accent",
      theme_success   = "success",
      theme_line      = "line"
    )

    ## Observe each colour picker and write to theme
    lapply(names(theme_inputs), function(input_id) {
      theme_key <- theme_inputs[[input_id]]
      shiny::observeEvent(input[[input_id]],
        {
          theme[[theme_key]] <- input[[input_id]]
          save_color_theme(shiny::reactiveValuesToList(theme), auth$user_dir)
        },
        ignoreInit = TRUE
      )
    })

    ## Palette dropdown
    shiny::observeEvent(input$theme_palette,
      {
        theme$palette <- input$theme_palette
        save_color_theme(shiny::reactiveValuesToList(theme), auth$user_dir)
      },
      ignoreInit = TRUE
    )

    ## Custom gradient colour pickers
    shiny::observeEvent(input$theme_palette_c1,
      {
        theme$palette_c1 <- input$theme_palette_c1
        save_color_theme(shiny::reactiveValuesToList(theme), auth$user_dir)
      },
      ignoreInit = TRUE
    )
    shiny::observeEvent(input$theme_palette_c2,
      {
        theme$palette_c2 <- input$theme_palette_c2
        save_color_theme(shiny::reactiveValuesToList(theme), auth$user_dir)
      },
      ignoreInit = TRUE
    )
    shiny::observeEvent(input$theme_palette_c3,
      {
        theme$palette_c3 <- input$theme_palette_c3
        save_color_theme(shiny::reactiveValuesToList(theme), auth$user_dir)
      },
      ignoreInit = TRUE
    )

    ## Reset button
    shiny::observeEvent(input$theme_reset, {
      defaults <- COLOR_THEME_DEFAULTS
      for (key in names(defaults)) {
        theme[[key]] <- defaults[[key]]
      }
      ## Also update the UI pickers to reflect defaults
      colourpicker::updateColourInput(session, "theme_primary", value = defaults$primary)
      colourpicker::updateColourInput(session, "theme_secondary", value = defaults$secondary)
      colourpicker::updateColourInput(session, "theme_neutral", value = defaults$neutral)
      colourpicker::updateColourInput(session, "theme_ns_color", value = defaults$ns_color)
      colourpicker::updateColourInput(session, "theme_bar_color", value = defaults$bar_color)
      colourpicker::updateColourInput(session, "theme_accent", value = defaults$accent)
      colourpicker::updateColourInput(session, "theme_success", value = defaults$success)
      colourpicker::updateColourInput(session, "theme_line", value = defaults$line)
      shiny::updateSelectInput(session, "theme_palette", selected = defaults$palette)
      colourpicker::updateColourInput(session, "theme_palette_c1", value = defaults$palette_c1)
      colourpicker::updateColourInput(session, "theme_palette_c2", value = defaults$palette_c2)
      colourpicker::updateColourInput(session, "theme_palette_c3", value = defaults$palette_c3)
      save_color_theme(shiny::reactiveValuesToList(theme), auth$user_dir)
    })

    ##----------------------------------------------------------------------
    ## new features
    ##----------------------------------------------------------------------

    newfeatures.RENDER <- reactive({
      newfeat <- markdown::markdownToHTML(
        file = file.path(OPG, "FEATURES.md"),
        fragment.only = TRUE
      )
      HTML(newfeat)
    })

    PlotModuleServer(
      "newfeatures",
      plotlib = "generic",
      func = newfeatures.RENDER,
      renderFunc = renderUI
    )

    ##----------------------------------------------------------------------
    ## packages
    ##----------------------------------------------------------------------
    
    packages.RENDER <- function() {
      pkg <- read.table(file.path(OPG, "RPackageLicenses.txt"), sep = "\t", header = TRUE)
      DT::datatable(pkg,
        rownames = FALSE,
        plugins = "scrollResize",
        fillContainer = TRUE,
        options = list(
          dom = "t",
          pageLength = 999,
          scrollResize = TRUE
        ),
        class = "compact hover"
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    TableModuleServer(
      "packages",
      func = packages.RENDER
    )

    ## Expose the runtime "Enable AI" switch to the app-level server, which
    ## show/hides the AI Studio + Copilot tabs in the root "app-sidebar"
    ## navset (see server.R). The deployment licence gate (opt$ENABLE_AI) is
    ## combined there.
    list(
      # Raw switch value (NULL before init) so the app-level server can treat
      # "not yet set" as on (default) and avoid a startup hide/show flash.
      enable_ai = shiny::reactive(input$enable_ai)
    )
  })
}
