##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ui.startupModal <- function(id, messages, title = NULL) {
  if (length(messages) == 0) {
    return(NULL)
  }

  header <- sapply(strsplit(messages, split = ":::"), function(m) m[[1]])
  messages2 <- sapply(messages, function(s) sub(".*:::", "", s))

  carousel_items <- list()
  for (i in 1:length(messages)) {
    tag1 <- bsutils::carouselItem(
      div(
        style = "height: 380px; margin-top: 0px;",
        class = "d-flex align-items-center justify-content-center",
        HTML(paste0(
          "<div><h4 class='modal-title text-center'>",
          header[[i]], "</h4>", messages2[[i]], "</div>"
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
  modal <- div(id = id, modal)
  return(modal)
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
