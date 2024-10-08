##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

WelcomeBoardInputs <- function(id) {
  return(NULL)
}


WelcomeBoardUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  pages <- list(
    '<h1 class="d-block w-100 text-center align-middle" style="height:100vh;line-height:90vh;">HELLO...</h1>',
    '<h1 class="d-block w-100 text-center align-middle" style="height:100vh;line-height:90vh;">WORLD!</h1>'
  )

  ## Pages for the slider
  mission.page <-
    div(
      class = "row welcome-slide",
      div(
        class = "col-md-12 text-center",
        shiny::tags$b("Our mission"),
        shiny::p("We love Biology. We love Big Data. Our mission is to create smart tools and make advanced omics analysis accessible to everyone. We believe that we can better understand Biology through Big Data, to find new cures and to accelerate the transition to data-driven precision medicine. Let’s together endeavour a world without cancer and complex diseases.")
      )
    )

  motto.page <-
    div(
      class = "row welcome-slide",
      div(
        class = "col-md-12 text-center",
        shiny::tags$b("Advanced omics analysis for everyone."), br(),
        "At BigOmics, we are focused on one thing — empowering biologists to easily visualize and understand their omics data. With Omics Playground you can analyze your omics data faster, better, easier with more fun. No coding required."
      )
    )

  credits.page <-
    div(
      class = "row welcome-slide",
      div(
        class = "col-md-12 text-center",
        shiny::tags$b("Proudly presented to you by"),
        shiny::p(
          "Ana Nufer, Antonino Zito, Axel Martinelli, Carson Sievert, Cédric Scherer, Gabriela Scorici, Griffin Seidel, Ivo Kwee, John Coene, Layal Abo Khayal, Marco Sciaini, Matt Leech, Mauro Miguel Masiero, Murat Akhmedov, Nick Cullen, Santiago Caño Muñiz, Shalini Pandurangan, Stefan Reifenberg, Xavier Escribà Montagut"
        )
      )
    )

  created.page <-
    div(
      class = "row welcome-slide",
      div(
        class = "col-md-12 text-center",
        shiny::tags$b("Created with love"), br(),
        "by BigOmics Analytics from Ticino, the sunny side of Switzerland.",
        br(), "© 2000-2024 BigOmics Analytics, Inc.", br(),
        shiny::a("www.bigomics.ch", href = "https://www.bigomics.ch")
      )
    )

  ## --------------------- page ------------------------------------------
  div(
    id = "welcome-page",
    ##    uiOutput(ns("notification")),
    div(
      class = "row",
      div(
        class = "col-md-12",
        br(),
        br(),
        div(shiny::textOutput(ns("welcome"), inline = TRUE), id = "welcome-text"),
        div(shiny::textOutput(ns("welcome2"), inline = TRUE), id = "welcome-subtext"),
        br(),
        br()
      )
    ),
    div(
      class = "row",
      style = "max-height:35vh;padding:0px 0px 20px 0px;vertical-align:bottom;",
      id = "welcome-buttons",
      div(
        class = "col-md-5",
        h3("I am new..."),
        shiny::actionButton(
          ns("btn_example_data"),
          label = "Load example dataset",
          class = "btn btn-outline-info welcome-btn"
        )
      ),
      div(
        class = "col-md-7",
        h3("I'm an existing user..."),
        shiny::actionButton(
          ns("btn_upload_new"),
          label = "Upload new data",
          class = "btn btn-outline-info welcome-btn"
        ),
        shiny::actionButton(
          ns("btn_load_data"),
          label = "Load from library",
          class = "btn btn-outline-primary welcome-btn"
        )
      )
    ),
    br(),
    br(),
    bs_carousel2(
      "welcome-carousel",
      wrap = TRUE, autostart = TRUE, fade = TRUE,
      interval = 10000,
      contents = list(mission.page, motto.page, created.page, credits.page)
    )
  )
}


WelcomeBoard <- function(id, auth, load_example, new_upload) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    bigdash.unloadSidebar()

    output$welcome <- shiny::renderText({
      shiny::req(auth$logged)
      if (!auth$logged) {
        return(NULL)
      }

      name <- auth$username
      if (is.null(name) || name %in% c("", NA)) name <- auth$email
      dbg("[WelcomeBoard] name = ", name)
      if (is.null(name) || name %in% c("", NA)) {
        welcome <- paste0("Welcome back...")
      } else {
        first.name <- getFirstName(name) ## in app/R/utils.R
        welcome <- paste0("Welcome back ", first.name, "...")
      }
      dbg("[WelcomeBoard] welcome = ", welcome)
      welcome
    })

    output$welcome2 <- shiny::renderText({
      shiny::req(auth$logged)
      if (!auth$logged) {
        return(NULL)
      }
      all.hello <- c(
        "Hello", "Salut", "Hola", "Pivet", "Ni hao", "Ciao", "Hi", "Hoi", "Hej",
        "Yassou", "Selam", "Hey", "Hei", "Grutzi", "Bonjour", "Jak się masz",
        "Namaste", "Salam", "Selamat", "Shalom", "Goeiedag", "Yaxshimusiz"
      )
      hello1 <- sample(all.hello, 1)
      paste0(hello1, "! What would you like to do today?")
    })

    observeEvent(input$btn_example_data, {
      if (is.null(load_example())) {
        load_example(1)
      } else {
        load_example(load_example() + 1)
      }
    })

    observeEvent(
      {
        input$btn_upload_new
      },
      {
        bigdash.selectTab(session, selected = "upload-tab")
      }
    )

    observeEvent(input$btn_load_data, {
      bigdash.openSettings(lock = TRUE)
      bigdash.openSidebar()
      bigdash.selectTab(session, "load-tab")
    })

    output$notification <- renderUI({
      div(bs_alert(HTML("Notification!"), style = "warning"), class = "p-5", style = "top: 10px !important; margin: 0 50px -40px 50px !important")
    })
  })
}
