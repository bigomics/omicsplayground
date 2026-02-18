##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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
    ## Plot Colors theme handling
    ## ----------------------------------------------------------------
    theme <- get_color_theme()

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
      shiny::observeEvent(input[[input_id]], {
        theme[[theme_key]] <- input[[input_id]]
      }, ignoreInit = TRUE)
    })

    ## Palette dropdown
    shiny::observeEvent(input$theme_palette, {
      theme$palette <- input$theme_palette
    }, ignoreInit = TRUE)

    ## Custom gradient colour pickers
    shiny::observeEvent(input$theme_palette_c1, {
      theme$palette_c1 <- input$theme_palette_c1
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$theme_palette_c2, {
      theme$palette_c2 <- input$theme_palette_c2
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$theme_palette_c3, {
      theme$palette_c3 <- input$theme_palette_c3
    }, ignoreInit = TRUE)

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
    })

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

    packages.RENDER <- function() {
      pkg <- read.table(file.path(OPG, "RPackageLicenses.txt"), sep = "\t", header = TRUE)
      DT::datatable(pkg,
        rownames = FALSE,
        plugins = "scrollResize",
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
  })
}
