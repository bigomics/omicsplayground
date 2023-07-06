##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

UserProfileBoard <- function(id, user) {
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

    observeEvent(user$logged(), {
      if (!user$logged()) {
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
      if (user$level() == "premium") {
        plan_class <- "success"
      }
      cl <- sprintf("badge badge-%s", plan_class)
      p(
        span("Subscription level", style = "color:grey;"),
        span(class = cl, tools::toTitleCase(user$level()))
      )
    })

    observeEvent(input$manage, {
      dbg("[UserBoard] !!! user$email() = ", user$email())
      dbg("[UserBoard] !!! user$stripe_id() = ", user$stripe_id())
      dbg("[UserBoard] !!! user$href = ", user$href())

      response <- httr::POST(
        "https://api.stripe.com/v1/billing_portal/sessions",
        body = list(
          customer = user$stripe_id(),
          return_url = user$href()
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
          Name   = user$name(),
          Email  = user$email(),
          Plan   = user$level(),
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
