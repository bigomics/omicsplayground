##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

WelcomeBoard <- function(id, auth, enable_upload, r_global) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    getFirstName <- reactive({
      name <- auth$name()
      if(is.null(name) || is.na(name) || name=='') name <- auth$email()
      if(is.null(name) || is.na(name) || name=='') {
        name <- "anonymous"
        name <- paste0("user",substring(session$token,1,4))
        return(name)
      }
      first.name <- strsplit(name, split = "[@ .]")[[1]][1]
      first.name <- paste0(
        toupper(substring(first.name, 1, 1)),
        substring(first.name, 2, nchar(first.name))
      )
      first.name
    })
    
    output$welcome <- shiny::renderText({
      name <- getFirstName()
      if(grepl("^user",name)) name <- ""
      if (name %in% c("", NA, NULL,"anonymous")) {
        welcome <- "Welcome back..."
      } else {
        welcome <- paste0("Welcome back ", name, "...")
      }
      welcome
    })

    observeEvent(input$btn_example_data, {
      r_global$load_example_trigger <- r_global$load_example_trigger + 1
    })

    observeEvent(input$btn_upload_data, {
      if (enable_upload) {
        bigdash.openSidebar()
        bigdash.selectTab(session, "upload-tab")
      } else {
        shinyalert::shinyalert(
          title = "Upload disabled",
          text = "Sorry, upload of new data is disabled for this account.",
          type = "warning",
          #
          closeOnClickOutside = FALSE
        )
      }
    })

    observeEvent(input$btn_load_data, {
      bigdash.openSettings(lock = TRUE)
      bigdash.openSidebar()
      bigdash.selectTab(session, "load-tab")
    })

    ## -------------------------------------------------------------
    ## Chatbox
    ## -------------------------------------------------------------
    if(opt$ENABLE_CHIRP) {
      shinyChatR::chat_server(
        "chatbox",
        db_file = file.path(OPG,"etc/chirp_data.db"),
        ##csv_file = file.path(OPG,"chirp_data.csv"),            
        chat_user = getFirstName
      )
    }
    
  })
}

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
        shiny::p("Ana Nufer, Axel Martinelli, Carson Sievert, Cédric Scherer, Gabriela Scorici, Ivo Kwee, John Coene, Layal Abo Khayal, Marco Sciaini, Matt Leech, Mauro Miguel Masiero, Murat Akhmedov, Nick Cullen, Stefan Reifenberg, Xavier Escribà Montagut")
      )
    )

  created.page <-
    div(
      class = "row welcome-slide",
      div(
        class = "col-md-12 text-center",
        shiny::tags$b("Created with love"), br(),
        "by BigOmics Analytics from Ticino, the sunny side of Switzerland.",
        br(), "© 2000-2023 BigOmics Analytics, Inc.", br(),
        shiny::a("www.bigomics.ch", href = "https://www.bigomics.ch")
      )
    )

  div.chirp <- NULL
  if(opt$ENABLE_CHIRP) {
    ## offcanvas chatbox 
    div.chirp <- bsutils::offcanvas(
      bsutils::offcanvasButton("Chirp!",id="chirp-button",
        style="width:auto;position:absolute;top:60px;right:30px;background-color:white;padding:5px 12px;font-size:14px;"),
      bsutils::offcanvasContent(
        .position = "end",
        bslib::card(
          full_screen = TRUE,
          style = "border-width: 0px;",
          height = "92vh",
          bslib::card_body(
            shinyChatR::chat_ui(ns("chatbox"),
              title = "Chirp with others on the Playground!",
              height='70vh', width='100%')
          )
        )
      )
    )
  }
  
  ## --------------------- page ------------------------------------------
  div(
    id = "welcome-page",
    div.chirp,
    div(
      class = "row",
      div(
        class = "col-md-12",
        br(),
        br(),
        div(shiny::textOutput(ns("welcome")), id = "welcome-text"),
        div("What would you like to do today?", id = "welcome-subtext"),
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
          ns("btn_upload_data"),
          label = "Upload new data",
          class = "btn btn-outline-info welcome-btn"
        ),
        shiny::actionButton(
          ns("btn_load_data"),
          label = "Use my saved data",
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
