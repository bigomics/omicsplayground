
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
