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

#' Modal with a tabset (navset) body.
#'
#' `tabs` is a named list, one element per tab: name = tab title, value =
#' the tab's UI content. A single-tab list skips the tab switcher entirely
#' and shows that one tab's content directly -- there's nothing to switch
#' between, and a lone tab header would just be visual noise.
#'
#' @param tabs Named list of tab title -> UI content.
#' @param header Modal header, e.g. a dataset name. `NULL` for none.
#' @param footer Modal footer. `NULL` for none.
#' @param size One of modalDialog2()'s sizes ("s", "m", "l", "xl",
#'   "fullscreen", "midscreen").
#' @param body_height CSS height (e.g. "65vh"), applied to *each tab's own
#'   content* (not the tab-header row above it) so every tab occupies
#'   exactly the same box and the modal's own size stays put as the user
#'   switches between them -- content scrolls (or, for an image built with
#'   `ui.imageTag(fit=TRUE, max_height=body_height)`, scales to fit) inside
#'   that box instead of growing/shrinking the dialog. Sized per tab content
#'   rather than around the whole navset so the tab-header row's own height
#'   doesn't eat into the budget an image is fit against -- pass the same
#'   value to both, or the image's cap and its actual box disagree and it
#'   scrolls instead of fitting. Only applied with more than one tab -- a
#'   single tab has nothing to switch between, so there's nothing to keep
#'   steady against. `NULL` to size each tab to its own content instead.
ui.showTabsetModal <- function(tabs, header = NULL, footer = NULL,
                               size = "l", body_height = "65vh",
                               easyClose = TRUE, fade = TRUE) {
  stopifnot(is.list(tabs), length(tabs) > 0)
  if (is.null(names(tabs)) || any(!nzchar(names(tabs)))) {
    stop("ui.showTabsetModal(): `tabs` must be a fully named list")
  }

  multi <- length(tabs) > 1
  if (multi && !is.null(body_height)) {
    tabs <- lapply(tabs, function(content) {
      div(style = paste0("height:", body_height, "; width:100%; overflow-y:auto;"), content)
    })
  }

  body <- if (!multi) {
    tabs[[1]]
  } else {
    do.call(bslib::navset_tab, lapply(names(tabs), function(nm) {
      bslib::nav_panel(title = nm, tabs[[nm]])
    }))
  }

  shiny::showModal(
    modalDialog2(
      body,
      header = header,
      footer = footer,
      size = size,
      easyClose = easyClose,
      fade = fade,
      ## Pin the dialog to its size class's max-width instead of letting it
      ## size to whichever tab's content happens to be showing (see
      ## modalDialog2()'s dialog_style docs).
      dialog_style = if (multi) "width:100%;" else NULL
    )
  )
}

#' @param dialog_style Extra CSS applied directly to the `.modal-dialog`
#'   element (in addition to the `size` class). The size classes only set a
#'   `max-width` -- `.modal-dialog` otherwise sizes to its content, which
#'   for a tabset can differ tab to tab (e.g. a narrower image vs. a
#'   full-width text panel) and visibly changes the dialog's width when
#'   switching. `dialog_style = "width:100%;"` (what `ui.showTabsetModal()`
#'   passes for a multi-tab modal) pins it to that max-width consistently,
#'   regardless of which tab is showing.
modalDialog2 <- function(
  ..., header = NULL, footer = modalButton("Dismiss"),
  size = c("m", "s", "l", "xl", "fullscreen", "midscreen"), easyClose = FALSE,
  fade = TRUE, dialog_style = NULL
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
      style = dialog_style,
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
