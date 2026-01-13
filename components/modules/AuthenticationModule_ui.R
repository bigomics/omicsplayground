##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## ================================================================================
## ================ AUTHENTICATION_MODULE UI FUNCTIONS ============================
## ================================================================================

splashLoginModal <- function(ns = NULL,
                             with.email = TRUE,
                             with.password = TRUE,
                             with.username = FALSE,
                             with.firebase = FALSE,
                             with.firebase_emailonly = FALSE,
                             with.link = FALSE,
                             link = NULL,
                             hide.password = TRUE,
                             id = "login",
                             button.text = "Log in",
                             cancel.text = "Cancel",
                             add.cancel = FALSE,
                             title = "Log in",
                             subtitle = "") {
  if (is.null(ns)) {
    ns2 <- function(e) {
      return(e)
    }
  }
  ns2 <- function(e) ns(paste0(id, "_", e))


  slogan <- list()
  slogan[[1]] <- c("Big Omics Data", "Isn't big anymore with BigOmics Playground")
  slogan[[2]] <- c("Great Discoveries", "Start on BigOmics Playground")
  slogan[[3]] <- c("Fasten Your Seat Belts!", "Hi-speed analytics")
  slogan[[4]] <- c("Do-it-yourself Omics Analytics", "Yes you can!")
  slogan[[5]] <- c("Twenty-Four Seven", "Your Playground doesn't go on coffee breaks")
  slogan[[6]] <- c("Analyze with confidence", "Be a data rockstar, a Freddie Mercury of omics!")
  slogan[[7]] <- c("Play-Explore-Discover", "Get deeper insights with BigOmics Playground")
  slogan[[8]] <- c("Skip the Queue", "Take the fast lane. Self-service analytics.")
  slogan[[9]] <- c("Look Ma! No help!", "I did it without a bioinformatician")
  slogan[[10]] <- c("Easy-peasy insight!", "Get insight from your data the easy way")
  slogan[[11]] <- c("Zoom-zoom-insight!", "Get faster insight from your data")
  slogan[[12]] <- c("Click-click-eureka!", "Owe yourself that <i>eureka!</i> moment")
  slogan[[13]] <- c("I Love Omics Data!", "Unleash your inner nerd with BigOmics Playground")
  slogan[[14]] <- c("More Omics Data", "Is all I want for my birthday")
  slogan[[15]] <- c("Keep Exploring", "Never stop discovering with BigOmics Playground")
  slogan[[16]] <- c("Real Bioinformaticians", "Do it with BigOmics Playground")
  slogan[[17]] <- c("Real Biologists", "Do it with BigOmics Playground")
  slogan[[18]] <- c("Ich bin doch nicht bl\u00F6d!", "Of course I use BigOmics Playground")
  slogan[[19]] <- c("Non sono mica scemo!", "Of course I use BigOmics Playground")
  slogan[[20]] <- c("The Unexplored Plan", "When you get into exploring, you realize that we live on a relatively unexplored plan. &ndash; E. O. Wilson")
  slogan[[21]] <- c("Explore More", "The more you explore, the more you learn and grow")
  slogan[[22]] <- c("Discover New Oceans", "Man cannot discover new oceans unless he has the courage to lose sight of the shore. &ndash; Andre Gide")
  slogan[[23]] <- c("Love Adventurous Life", "Be passionately curious about exploring new adventures. &ndash; Lailah Gifty Akita")
  slogan[[24]] <- c("Succes is Exploration", "Find the unknown. Learning is searching. Anything else is just waiting. &ndash; Dale Daute")
  slogan[[25]] <- c("Look Ma! No help!", "I did it without a bioinformagician")
  slogan[[26]] <- c("May the Force of Omics be with you", "Train hard youngling, one day a master you become")

  slogan <- sample(slogan, 1)[[1]]
  slogan.len <- nchar(paste(slogan, collapse = " "))
  if (slogan.len < 80) slogan[1] <- paste0("<br>", slogan[1])
  splash.title <- shiny::div(
    id = "splash-title",
    class = "text-white",
    shiny::div(shiny::HTML(slogan[1]), style = "font-size:3rem;font-weight:700;line-height:1em;width:130%;"),
    shiny::div(shiny::HTML(slogan[2]), style = "font-size:1.6rem;line-height:1.1em;margin-top:0.5em;width:130%;")
  )

  div.password <- div()
  div.email <- div()
  div.username <- div()
  div.firebase <- div()
  div.title <- div()
  div.subtitle <- div()
  div.link <- div()

  if (with.email) {
    div.email <- div(
      id = "splash-email",
      ## tags$input(
      ##   type = "email",
      ##   id = ns2("login_email"),
      ##   placeholder = "your email",
      ##   autocomplete = "email",
      ##   class = "form-control shiny-bound-input shinyjs-resettable"
      ## )
      textInput(ns2("email"), NULL, placeholder = "your email")
    )
  }
  if (with.username) {
    div.username <- div(
      id = "splash-username",
      textInput(ns2("username"), NULL, placeholder = "your username")
    )
  }
  if (with.password) {
    if (hide.password) {
      div.password <- div(
        id = "splash-password",
        passwordInput(ns2("password"), NULL, placeholder = "your password")
      )
    } else {
      div.password <- div(
        id = "splash-password",
        textInput(ns2("password"), NULL, placeholder = "your password")
      )
    }
  }
  if (with.firebase) {
    div.firebase <- div(
      class = "card",
      div(
        class = "card-body",
        h1(
          title,
          class = "card-title pb-2"
        ),
        div(
          subtitle
        ),
        textInput(
          ns2("emailInput"),
          "",
          placeholder = "Your email",
          width = "100%"
        ),
        actionButton(
          ns2("emailSubmit"),
          "Send link",
          class = "btn-warning"
        ),
        p(
          id = "emailFeedbackShow"
        ),
        hr(style = "color:#888;opacity:1;margin-top:30px;"),
        h5(
          "or",
          class = "text-center pb-4 pt-1",
          style = "margin-top:-35px;background:white;width:50px;margin-left:auto;margin-right:auto;"
        ),
        div(
          class = "social-button google-button",
          actionLink(
            ns2("launchGoogle"),
            HTML("&nbsp; Sign in with Google"),
            icon = shiny::icon("google", style = "font-size:18px;")
          )
        ),
        div(
          class = "social-button facebook-button",
          actionLink(
            ns2("launchFacebook"),
            HTML("&nbsp; Sign in with Facebook"),
            icon = shiny::icon("facebook", style = "font-size:18px;")
          )
        ),
        ## div(
        ##   class = "social-button apple-button",
        ##   actionLink(
        ##     ns2("launchApple"),
        ##     HTML("&nbsp; Sign in with Apple"),
        ##     icon = shiny::icon("apple", style="font-size:18px;")
        ##   )
        ## ),
        div(
          class = "social-button twitter-button",
          actionLink(
            ns2("launchTwitter"),
            HTML("&nbsp; Sign in with Twitter"),
            icon = shiny::icon("twitter", style = "font-size:18px;")
          )
        )
      )
    )
  }
  if (with.firebase_emailonly) {
    div.firebase <- div(
      class = "card",
      div(
        class = "card-body",
        h1(
          title,
          class = "card-title pb-2"
        ),
        div(
          subtitle
        ),
        textInput(
          ns2("emailInput"),
          "",
          placeholder = "Your email",
          width = "100%"
        ),
        actionButton(
          ns2("emailSubmit"),
          "Send link",
          class = "btn-warning"
        ),
        p(
          id = "emailFeedbackShow"
        )
      )
    )
  }

  if (!is.null(title) && title != "") {
    div.title <- div(
      id = ns2("splash-login-title"),
      class = "pb-3",
      h1(title, style = "color:black;line-height:1em;")
    )
  }

  if (!is.null(subtitle) && subtitle != "") {
    div.subtitle <- div(
      id = ns2("splash-login-subtitle"),
      class = "pt-0 pb-2",
      h6(subtitle, style = "color:black;font-weight:400;")
    )
  }

  if (add.cancel) {
    div.button <- div(
      id = ns2("splash-buttons"),
      class = "pt-2",
      shiny::fillRow(
        flex = c(1, NA, NA, 1),
        br(),
        actionButton(ns2("cancel_btn"), cancel.text,
          class = "btn-light btn-lg",
          style = "margin: 4px;"
        ),
        actionButton(ns2("submit_btn"), button.text,
          class = "btn-warning btn-lg",
          style = "margin: 4px;"
        ),
        br()
      )
    )
  } else {
    div.button <- div(
      id = "splash-buttons",
      class = "pt-2",
      actionButton(ns2("submit_btn"), button.text, class = "btn-warning btn-xl")
    )
  }

  if (with.link) {
    div.link <- div(
      id = ns2("link"),
      class = "pt-2",
      tags$a(
        id = ns2("link_btn"),
        href = link,
        target = "_blank",
        button.text,
        class = "btn btn-warning btn-xl",
        role = "button"
      )
    )
    div.button <- div()
  }

  ## splash.panel=div();ns=function(x)
  splash.content <- NULL
  if (with.firebase || with.firebase_emailonly) {
    splash.content <- div.firebase
  } else {
    splash.content <- shiny::wellPanel(
      ## style = "padding: 40px 20px; background-color: #ffffff22;",
      style = "padding: 35px 25px; background-color:white; color:black;",
      id = ns2("splash-login"),
      div.title,
      div.subtitle,
      div.username,
      div.email,
      div.password,
      div.button,
      div.link
    )
  }

  body <- div(
    id = ns2("splash-content"),
    splash.content
  )

  ## display the number of active sessions
  num_sessions <- paste0(length(ACTIVE_SESSIONS), "/", MAX_SESSIONS)

  footer <- div(
    id = ns2("splash-footer"),
    style = "position: absolute; bottom:5px; left:10px; color:#ffffff88; font-size:0.85em;",
    getAppVersion(add.auth = TRUE),
    div(num_sessions, style = "padding-left:1em; display: inline;")
  )

  m <- splashScreen(title = splash.title, body = body, ns = ns2, footer2 = footer)
  return(m)
}

splashscreen.buttons <- function(ns) {
  tagList(
    shiny::tags$a(
      shiny::img(
        id = ns("splash-logo2"),
        class = "splash-logo2",
        src = "static/bigomics-logo.png"
      ),
      href = "https://www.bigomics.ch",
      target = "_blank"
    ),
    div(
      class = "btn-group",
      role = "group",
      div(
        class = "btn-group",
        role = "group",
        tags$button(
          "Get support",
          id = "splash-toggle-support",
          type = "button",
          class = "btn btn-outline-primary dropdown-toggle",
          `data-bs-toggle` = "dropdown",
          `aria-expanded` = "false"
        ),
        tags$ul(
          class = "dropdown-menu",
          `aria-labelledby` = "splash-toggle-support",
          tags$li(
            shiny::tags$a(
              "Watch tutorials",
              class = "dropdown-item",
              href = "https://www.youtube.com/channel/UChGASaLbr63pxmDOeXTQu_A",
              target = "_blank"
            )
          ),
          tags$li(
            shiny::tags$a(
              "Read documentation",
              class = "dropdown-item",
              href = "https://omicsplayground.readthedocs.io",
              target = "_blank"
            ),
          ),
          tags$li(
            shiny::tags$a(
              "User forum",
              class = "dropdown-item",
              href = "https://groups.google.com/d/forum/omicsplayground",
              target = "_blank"
            )
          )
        )
      ),
      div(
        class = "btn-group",
        role = "group",
        tags$button(
          "I'm a developer",
          id = ns("splash-toggle-dev"),
          type = "button",
          class = "btn btn-outline-primary dropdown-toggle",
          `data-bs-toggle` = "dropdown",
          `aria-expanded` = "false"
        ),
        tags$ul(
          class = "dropdown-menu",
          `aria-labelledby` = "splash-toggle-dev",
          tags$li(
            shiny::tags$a(
              "Get the source",
              class = "dropdown-item",
              href = "https://github.com/bigomics/omicsplayground",
              target = "_blank"
            )
          ),
          tags$li(
            shiny::tags$a(
              "Docker image",
              class = "dropdown-item",
              href = "https://hub.docker.com/r/bigomics/omicsplayground",
              target = "_blank"
            )
          ),
          tags$li(
            shiny::tags$a(
              "Buy us a coffee!",
              class = "dropdown-item",
              href = "https://www.buymeacoffee.com/bigomics",
              target = "_blank"
            )
          )
        )
      )
    )
  )
}

splashScreen <- function(title, body, ns = NULL, easyClose = FALSE, fade = FALSE,
                         buttons = TRUE, footer = NULL, footer2 = NULL) {
  if (is.null(ns)) {
    ns <- function(e) {
      return(e)
    }
  }

  div.buttons <- shiny::modalButton("Dismiss")
  if (buttons) {
    div.buttons <- splashscreen.buttons(ns = ns)
  }

  ## return modalDialog
  m <- modalDialog2(
    id = ns("splash-fullscreen"),
    class = "bg-primary",
    header = div.buttons,
    shiny::div(
      class = "row",
      shiny::div(
        class = "col-md-4 offset-md-2",
        title,
        br(),
        br(),
        shiny::img(src = "static/mascotte-sc.png", class = "img-fluid", id = "splash-image"),
      ),
      shiny::div(
        class = "col-md-3 offset-md-2",
        shiny::div(
          id = "splash-panel",
          body,
          div(textOutput(ns("warning")), style = "color:white;font-size:1.2em;padding-top:8px;line-height:1.1em;"),
        ),
      )
    ),
    footer2,
    footer = footer,
    size = "fullscreen",
    easyClose = easyClose,
    fade = fade
  ) ## end of modalDialog

  return(m)
}


## ================================================================================
## ================================= END OF FILE ==================================
## ================================================================================
