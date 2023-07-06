##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

UserSettingsBoard <- function(id, user) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    dbg("[UserSettingsBoard] >>> initializing User Settings...")

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
