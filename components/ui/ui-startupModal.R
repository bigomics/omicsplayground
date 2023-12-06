##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ui.startupModal <- function(id, messages) {
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

  modal <- bsutils::modal(
    id = id,
    size = "lg",
    bsutils::modalHeader(
      bsutils::modalTitle("What's new on the Playground"),
      style = "background-color: #b6d3E888"
    ),
    bsutils::modalBody(
      shiny::tagAppendAttributes(
        do.call(
          function(...) {
            bsutils::carousel(
              ...,
              indicators = TRUE, controls = TRUE
            )
          },
          carousel_items
        ),
        `data-bs-interval` = "12000"
      ),
      style = "background-color: #b6d3E822"
    ),
    centered = TRUE
  )

  modal <- shiny::tagAppendAttributes(
    modal,
    `data-bs-backdrop` = "false"
  )

  return(modal)
}
