##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

ui.showSmallModal <- function(msg = "Please wait...", timer = 0) {
  shiny::showModal(shiny::modalDialog(
    title = NULL,
    shiny::HTML("<br><center><p>", msg, "</p></center>"),
    footer = NULL,
    size = "s", easyClose = TRUE, fade = FALSE
  ))
  if (timer > 0) {
    shinyjs::delay(timer, shiny::removeModal())
  }
}


pgx.showSmallModal <- function(msg = "Please wait...") {
  shiny::showModal(shiny::modalDialog(
    title = NULL,
    shiny::HTML("<br><center><p>", msg, "</p></center>"),
    footer = NULL,
    size = "s",
    easyClose = FALSE,
    fade = FALSE
  ))
}

pgx.showSmallModal2 <- function(msg = "Please wait...", easyClose = TRUE,
                                footer = modalButton("Dismiss")) {
  shiny::showModal(shiny::modalDialog(
    title = NULL,
    shiny::HTML("<br><center><p>", msg, "</p></center>"),
    footer = div(footer, class = "text-center"),
    size = "s",
    easyClose = easyClose,
    fade = FALSE
  ))
}

modalDialog2 <- function(
  ..., header = NULL, footer = modalButton("Dismiss"),
  size = c("m", "s", "l", "xl", "fullscreen", "midscreen"), easyClose = FALSE, fade = TRUE
) {
  size <- match.arg(size)
  backdrop <- if (!easyClose) {
    "static"
  }
  keyboard <- if (!easyClose) {
    "false"
  }
  div(
    id = "shiny-modal", class = "modal", class = if (fade) {
      "fade"
    }, tabindex = "-1", `data-backdrop` = backdrop,
    `data-bs-backdrop` = backdrop, `data-keyboard` = keyboard,
    `data-bs-keyboard` = keyboard, div(
      class = "modal-dialog",
      class = switch(size,
        s = "modal-sm",
        m = NULL,
        l = "modal-lg",
        xl = "modal-xl",
        fullscreen = "modal-fullscreen",
        midscreen = "modal-midscreen"
      ),
      div(
        class = "modal-content",
        if (!is.null(header)) {
          div(class = "modal-header", header)
        }, div(class = "modal-body", ...),
        if (!is.null(footer)) {
          div(class = "modal-footer", footer)
        }
      )
    ), tags$script(HTML("if (window.bootstrap && !window.bootstrap.Modal.VERSION.match(/^4\\./)) {\n         var modal = new bootstrap.Modal(document.getElementById('shiny-modal'));\n         modal.show();\n      } else {\n         $('#shiny-modal').modal().focus();\n      }"))
  )
}
