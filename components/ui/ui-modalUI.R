##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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

modalTrigger <- function(
    id,
    target,
    ...,
    class = "") {
  class <- sprintf(
    "btn %s",
    class
  )

  tags$a(
    id = id,
    ...,
    `data-bs-toggle` = "modal",
    `data-bs-target` = sprintf("#%s", target),
    class = class
  )
}

modalUI <- function(
    id,
    title,
    ...,
    size = c("default", "sm", "lg", "xl", "fullscreen"),
    track_open = FALSE,
    footer = tags$div(
      class = "modal-footer",
      tags$button(
        type = "button",
        class = "btn btn-secondary",
        `data-bs-dismiss` = "modal",
        "Close"
      )
    )) {
  size <- match.arg(size)

  size_cl <- switch(size,
    "sm" = "modal-sm",
    "lg" = "modal-lg",
    "xl" = "modal-xl",
    "fullscreen" = "modal-fullscreen",
    ""
  )

  modal <- tags$div(
    class = "modal fade",
    id = id,
    tabindex = "-1",
    `aria-labelledby` = "exampleModalLabel",
    `aria-hidden` = "true",
    tags$div(
      class = sprintf("modal-dialog %s", size_cl),
      tags$div(
        class = "modal-content",
        tags$div(
          class = "modal-header",
          tags$div(
            class = "modal-title",
            title
          ),
          tags$button(
            type = "button",
            class = "btn-close",
            `data-bs-dismiss` = "modal",
            `aria-label` = "Close"
          )
        ),
        tags$div(
          class = "modal-body",
          ...
        ),
        footer
      )
    )
  )

  if (track_open) {
    tagList(
      modal,
      tags$script(sprintf(
        "$('#%s').on('shown.bs.modal hidden.bs.modal', function(e) {
          console.log('Modal event:', e.type, 'for modal:', '%s');
          Shiny.setInputValue('%s_is_open', e.type === 'shown');
        });",
        id,
        id,
        id
      ))
    )
  } else {
    modal
  }
}

modalDialog2 <- function(
    ..., header = NULL, footer = modalButton("Dismiss"),
    size = c("m", "s", "l", "xl", "fullscreen"), easyClose = FALSE, fade = TRUE) {
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
        fullscreen = "modal-fullscreen"
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
