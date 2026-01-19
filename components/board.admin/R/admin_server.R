##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

AdminPanelBoard <- function(id, auth) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    dbg("[AdminPanelBoard] >>> initializing AdminBoard...")

    ## Check if user is admin - this is a critical security check
    is_admin <- reactive({
      isTRUE(auth$ADMIN)
    })

    shiny::observeEvent(input$board_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Admin Panel</strong>"),
        shiny::HTML(
          "The Admin Panel provides administrative functions for managing
          the OmicsPlayground platform. This panel is only accessible to
          users with admin privileges."
        ),
        easyClose = TRUE, size = "l"
      ))
    })

    ## Admin info display
    output$admin_info <- renderUI({
      req(is_admin())
      tagList(
        p(strong("Current Admin: "), auth$username),
        p(strong("Email: "), auth$email),
        p(strong("Admin Status: "), span(class = "badge bg-success", "Active"))
      )
    })

    ## System status display
    output$system_status <- renderUI({
      req(is_admin())
      
      ## Get some basic system info
      r_version <- paste0(R.version$major, ".", R.version$minor)
      platform <- R.version$platform
      
      tagList(
        p(strong("R Version: "), r_version),
        p(strong("Platform: "), platform),
        p(strong("Working Directory: "), getwd())
      )
    })

    ## User statistics table
    output$user_stats <- renderTable(
      {
        req(is_admin())
        dbg("[AdminBoard::user_stats] renderTable")
        
        ## Get user directories as a proxy for registered users
        user_dirs <- list.dirs(PGX.DIR, full.names = FALSE, recursive = FALSE)
        user_dirs <- grep("@", user_dirs, value = TRUE)
        
        data.frame(
          Metric = c("Total User Directories", "Data Directory"),
          Value = c(length(user_dirs), PGX.DIR),
          check.names = FALSE
        )
      },
      width = "100%",
      striped = TRUE
    )


  })
}
