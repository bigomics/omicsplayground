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

    newfeatures.RENDER <- reactive({
      feature_file <- file.path(OPG, "FEATURES.md")
      newfeat <- xfun::read_utf8(feature_file)
      HTML(opg_markdown_to_html(paste(newfeat, collapse = "\n")))
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
