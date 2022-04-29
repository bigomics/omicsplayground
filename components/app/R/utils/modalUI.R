
modalTrigger <- function(
  id,
  target,
  ...,
  class = ""
) {
  class <- sprintf(
    "btn %s",
    class
  )

  tags$button(
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
  size = c("default", "sm", "lg", "xl","fullscreen"),
  footer = tags$div(
    class = "modal-footer",
    tags$button(
      type = "button",
      class = "btn btn-secondary",
      `data-bs-dismiss` = "modal",
      "Close"
    )
  )
){
  size <- match.arg(size)

  size_cl <- switch(
      size,
      "sm" = "modal-sm",
      "lg" = "modal-lg",
      "xl" = "modal-xl",
      "fullscreen" = "modal-fullscreen",
      ""
  )
  
  tags$div(
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
          tags$h5(
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
}


modalDialog2 <- function (..., title = NULL, footer = modalButton("Dismiss"), 
    size = c("m", "s", "l", "xl","fullscreen"), easyClose = FALSE, fade = TRUE) 
{
    size <- match.arg(size)
    backdrop <- if (!easyClose) 
        "static"
    keyboard <- if (!easyClose) 
        "false"
    div(id = "shiny-modal", class = "modal", class = if (fade) 
        "fade", tabindex = "-1", `data-backdrop` = backdrop, 
        `data-bs-backdrop` = backdrop, `data-keyboard` = keyboard, 
        `data-bs-keyboard` = keyboard, div(class = "modal-dialog", 
            class = switch(size, s = "modal-sm", m = NULL, l = "modal-lg", 
                           xl = "modal-xl", fullscreen = "modal-fullscreen"),
            div(class = "modal-content", 
                if (!is.null(title)) 
                  div(class = "modal-header", tags$h4(class = "modal-title", 
                    title)), div(class = "modal-body", ...), 
                if (!is.null(footer)) 
                  div(class = "modal-footer", footer))), tags$script(HTML("if (window.bootstrap && !window.bootstrap.Modal.VERSION.match(/^4\\./)) {\n         var modal = new bootstrap.Modal(document.getElementById('shiny-modal'));\n         modal.show();\n      } else {\n         $('#shiny-modal').modal().focus();\n      }")))
}
