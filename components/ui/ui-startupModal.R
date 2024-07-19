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
        style = "height: 360px;",
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
      style = "background-color: #f0f9fd; margin-bottom: 0px;"
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
    #    tags$style(".modal-dialog {width: 720px;}"),
    #    tags$style(".modal-content {background-color: #f0f9fd;}"),
    #    tags$style(".modal-header {padding: 0px;}")
  )
  modal <- div(id = id, modal)
  return(modal)
}
