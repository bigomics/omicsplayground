##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

UserSettingsBoard <- function(id, auth, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    dbg("[UserSettingsBoard] >>> initializing User Settings...")

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

#    output$newfeatures <- renderUI({
#      newfeat <- markdown::markdownToHTML(file = file.path(OPG, "FEATURES.md"),
#                                          fragment.only = TRUE)
#      HTML(newfeat)
#    })

    newfeatures.RENDER <- reactive({
        newfeat <- markdown::markdownToHTML(file = file.path(OPG, "FEATURES.md"),
                                            fragment.only = TRUE)
        HTML(newfeat)
    })

    PlotModuleServer(
      "newfeatures",
      plotlib = "generic",
      func = newfeatures.RENDER,
      renderFunc = renderUI
    )
    
    packages.RENDER <- function() {
      pkg <- read.table(file.path(OPG, "RPackageLicenses.txt"),sep="!",header=TRUE)[,1]
      pkg <- strsplit(gsub("[ ]{2,}",",", trimws(pkg)),split=",")
      pkg <- do.call(rbind, pkg)
      colnames(pkg) <- c("Package","Version","License")
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
    
    res <- list(
      enable_beta = reactive({
        as.logical(input$enable_beta)
      }),
      enable_info = reactive({
        as.logical(input$enable_info)
      })
    )
    return(res)
  })
}
