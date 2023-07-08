##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

UserProfileBoard <- function(id, auth) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    dbg("[UserProfileBoard] >>> initializing UserBoard...")

    shiny::observeEvent(input$board_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>User Profile</strong>"),
        shiny::HTML(
          "The User Profile page allows you to view and change your
                subscription plan and to view the latest news about application
                development."
        ),
        easyClose = TRUE, size = "l"
      ))
    })

    observeEvent(auth$logged, {
      if (!auth$logged) {
        return()
      }

      removeModal()

      session$sendCustomMessage(
        "get-subs",
        list(
          ns = ns(NULL)
        )
      )
    })

    output$plan <- renderUI({
      plan_class <- "info"
      if (auth$level == "premium") {
        plan_class <- "success"
      }
      cl <- sprintf("badge badge-%s", plan_class)
      p(
        span("Subscription level", style = "color:grey;"),
        span(class = cl, tools::toTitleCase(auth$level))
      )
    })

    observeEvent(input$manage, {
      dbg("[UserBoard] !!! auth$email = ", auth$email)
      dbg("[UserBoard] !!! auth$stripe_id = ", auth$stripe_id)
      dbg("[UserBoard] !!! auth$href = ", auth$href)

      response <- httr::POST(
        "https://api.stripe.com/v1/billing_portal/sessions",
        body = list(
          customer = auth$stripe_id,
          return_url = auth$href
        ),
        httr::authenticate(
          Sys.getenv("OMICS_STRIPE_KEY"),
          ""
        ),
        encode = "form"
      )

      httr::warn_for_status(response)
      content <- httr::content(response)
      session$sendCustomMessage("manage-sub", content$url)
    })

    output$userdata <- renderTable(
      {
        dbg("[UserBoard::userdata]  renderDataTable")
        cl <- "badge badge-info"
        values <- c(
          Name   = auth$username,
          Email  = auth$email,
          Plan   = auth$level,
          Start  = "",
          End    = "",
          Status = "active"
        )
        values[which(values == "")] <- "(not set)"
        data.frame(" " = names(values), "  " = values, check.names = FALSE)
      },
      width = "400px",
      striped = TRUE
    )

    output$news <- renderUI({
      news <- markdown::markdownToHTML(file = file.path(OPG, "VERSION"), fragment.only = TRUE)
      HTML(news)
    })
  })
}
