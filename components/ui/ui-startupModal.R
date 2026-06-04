##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ui.startupModal <- function(...) {
  ui.showStartupModal(...)
}

ui.showStartupModal <- function(title = "BigOmics Highlights", ...) {  

  ## read startup messages
  msg_file <- file.path(ETC, "MESSAGES")
  if (!file.exists(msg_file)) return(NULL)

  msg <- readLines(msg_file)
  msg <- msg[msg != "" & substr(msg, 1, 1) != "#"]
  #msg <- c(msg[[1]], sample(msg, min(4, length(msg))))
  msg <- sample(msg, min(4, length(msg)))
  
  header <- sapply(strsplit(msg, split = ":::"), function(m) m[[1]])
  msg2 <- sapply(msg, function(s) sub(".*:::", "", s))

  carousel_items <- list()
  for (i in 1:length(msg)) {
    tag1 <- bsutils::carouselItem(
      div(
        style = "height: 380px; margin-top: 0px;",
        class = "d-flex align-items-center justify-content-center",
        HTML(paste0(
          "<div><h4 class='modal-title text-center'>",
          header[[i]], "</h4>", msg2[[i]], "</div>"
        ))
      ),
      class = "p-2"
    )
    carousel_items[[i]] <- tag1
  }

  modal <- shiny::modalDialog(
    size = "l",
    title = NULL,
    footer = NULL,
    bsutils::modalHeader(
      bsutils::modalTitle(""),
      style = "background-color: #f0f9fd; margin-bottom: 0px; padding: 0px;"
    ),
    do.call(
      function(...) {
        bsutils::carousel(
          ...,
          id = "opg-welcome-carousel",
          indicators = TRUE,
          controls = TRUE
        )
      },
      carousel_items
    ),
    easyClose = TRUE
  )
  shiny::showModal(
    div(id = "startup_modal", modal)
  )
}

ui.showCartoonModal <- function(msg = "Loading data...", img.path = "www/cartoons") {
  cartoon_list <- list(
    list(slogan = "Visual analytics. See and understand", img = "data-graph-wisdom.jpg"),
    list(slogan = "Fasten your seat belts. Accelerated discovery", img = "cartoon-speedup.jpg"),
    list(slogan = "Analytics anywhere. Anytime.", img = "cartoon-cloudservice.jpg"),
    list(slogan = "Analyze with confidence. Be a rockstar", img = "bigomics-rockstar3.jpg"),
    list(slogan = "Fast track your Bioinformatics", img = "selfservice-checkout2.png"),
    list(slogan = "Integrate more. Dig deeper", img = "cartoon-integration.jpg"),
    list(slogan = "Your analysis doesn't take coffee breaks", img = "gone-for-coffee.png"),
    list(slogan = "Too much data? Help yourself", img = "cartoon-datahelp2.jpg"),
    list(slogan = "Big Friendly Omics", img = "big-friendly-omics1.jpg"),
    list(slogan = "Big Data meets Biology", img = "bigdata-meets.png")
  )

  randomCartoon <- function() {
    cartoon <- sample(cartoon_list, 1)[[1]]
    cartoon$img2 <- paste0("static/cartoons/", cartoon$img)
    cartoon$img <- file.path(img.path, cartoon$img)
    cartoon
  }

  toon <- randomCartoon()
  shiny::showModal(shiny::modalDialog(
    title = shiny::div(shiny::h2(toon$slogan), shiny::p("with Omics Playground"), style = "text-align:center;"),
    shiny::img(src = toon$img2, class = "img-fluid"),
    footer = fillRow(flex = c(1, NA, 1), " ", msg, " "),
    size = "l",
    easyClose = FALSE,
    fade = TRUE
  ))
}

ui.showImageModal <- function(img, title, footer='', width=1088) {
  imgfile = tempfile(fileext=".png")
  png::writePNG(img, target = imgfile)
  hratio <- dim(img)[1] / dim(img)[2]
  header <- NULL
  if(!is.null(title) && title != "") header <- shiny::h3(HTML(title))
  shiny::showModal(
    modalDialog2(
      shiny::div(shiny::img(src = base64enc::dataURI(file=imgfile),
        width = width, height = hratio*width)),
      header = header,
      footer = div(HTML(footer),
        style='display:inline; font-size:1.2em; line-height:1.1em'),
      size = "xl",
      easyClose = TRUE,
      fade = TRUE
    )
  )
}

ui.showAboutModal <- function() {
  authors <- c(
    "Ana Nufer, Antonino Zito, Axel Martinelli, Carson Sievert, Cédric Scherer, Gabriela Scorici, Griffin Seidel, Ivo Kwee, John Coene, Jonathan Manson-Hennig, Layal Abo Khayal, Marco Sciaini, Matt Leech, Mauro Miguel Masiero, Murat Akhmedov, Nick Cullen, Santiago Caño Muñiz, Shalini Pandurangan, Stefan Reifenberg, Xavier Escribà Montagut"
  )
  authors <- paste(sort(authors), collapse = ", ")
  
  shiny::showModal(
    shiny::modalDialog(
      div(
        h2("Omics Playground"),
        h5(VERSION),
        h5("Advanced omics analysis for everyone"), br(), br(),
        p("Created with love and proudly presented to you by BigOmics Analytics from Ticino, the sunny side of Switzerland."),
        p(tags$a(href="https://www.bigomics.ch", "www.bigomics.ch", target="_blank")),
        style = "text-align:center; line-height: 1em;"
      ),
      footer = div(
        "© 2000-2026 BigOmics Analytics, Inc.",
        br(), br(),
        paste("Credits:", authors),
        style = "font-size: 0.8em; line-height: 0.9em; text-align:center;"
      ),
      size = "m",
      easyClose = TRUE,
      fade = FALSE
    )
  )
}
