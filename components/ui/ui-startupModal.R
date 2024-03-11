##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ui.startupModal <- function(id, messages, title) {
  if (length(messages) == 0) {
    return(NULL)
  }

  carousel_items <- list()
  for (i in 1:length(messages)) {
    tag1 <- bsutils::carouselItem(
      div(
        style = "height: 360px;",
        class = "d-flex align-items-center justify-content-center",
        HTML(paste0("<div>", messages[[i]], "</div>"))
      ),
      class = "p-4"
    )
    carousel_items[[i]] <- tag1
  }

    modal <- shiny::modalDialog(
      size = "l",
      title = NULL,
      footer = NULL,
      bsutils::modalHeader(
        bsutils::modalTitle(title),
        style = "background-color: #b6d3E888"),
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
      easyClose = TRUE,
      # tag$style can be used to change background color
      tags$style(".modal-content {background-color: #b6d3E888;}")
    )
    return(modal)
  }
