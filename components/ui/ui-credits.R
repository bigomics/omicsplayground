##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

creditsCarousel <- function() {

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
          "Ana Nufer, Antonino Zito, Axel Martinelli, Carson Sievert, Cédric Scherer, Gabriela Scorici, Griffin Seidel, Ivo Kwee, John Coene, Jonathan Manson-Hennig, Layal Abo Khayal, Marco Sciaini, Matt Leech, Mauro Miguel Masiero, Murat Akhmedov, Nick Cullen, Santiago Caño Muñiz, Shalini Pandurangan, Stefan Reifenberg, Xavier Escribà Montagut"
        )
      )
    )

  created.page <-
    div(
      class = "row welcome-slide",
      div(
        class = "col-md-12 text-center",
        shiny::tags$b("Created with love by"), br(),
        "BigOmics Analytics from Ticino, the sunny side of Switzerland.",
        br(), "© 2000-2025 BigOmics Analytics, Inc.", br(),
        shiny::a("www.bigomics.ch", href = "https://www.bigomics.ch")
      )
    )
  
  bs_carousel2(
    "welcome-carousel",
    wrap = TRUE, autostart = TRUE, fade = TRUE,
    interval = 10000,
    contents = list(mission.page, motto.page, created.page, credits.page)
  )

}
