##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Various Bootstrap goodies.
##
##


bs_alert <- function(m="alert!") {
  shiny::tags$div(
    class = "alert alert-primary alert-dismissible fade show",
    role = "alert",
    m,
    shiny::tags$button(
      type = "button",
      class = "btn-close",
      `data-bs-dismiss` = "alert",
      `aria-label` = "Close",
      shiny::tags$span(
        `aria-hidden` = "true"
      )
    )
  )
}
