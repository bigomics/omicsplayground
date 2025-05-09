##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Various Bootstrap goodies.
##
##

bs_alert <- function(text, conditional = TRUE, style = "primary", closable = TRUE,
                     translate = TRUE, translate_js = TRUE, html = TRUE) {
  id <- bigdash:::make_id()
  if (conditional) {
    btn.class <- "btn-close btn-close-bs-conditional"
  } else {
    btn.class <- "btn-close"
  }

  ##  text <- list(...)

  if (translate) {
    text <- tspan(text, js = translate_js)
  }

  if (html) {
    text <- shiny::HTML(text)
  }

  shiny::tags$div(
    id = id,
    class = paste0("alert alert-", style, " alert-dismissible fade show"),
    role = "alert",
    text,
    if (closable) {
      shiny::tags$button(
        # Use display: none; instead of official boostrap close button.
        # If not, the element interfers with the bslib::layout_columns
        # and we get extra gap on top. The second part of the `onclick` is
        # to also close the box (hide it) when it is not placed inside a
        # bslib::layout_columns
        onclick = paste0('$("#', id, ' button").closest(".bslib-gap-spacing.html-fill-container").css("display", "none");$("#', id, ' button").parent().css("display", "none");'),
        type = "button",
        class = btn.class,
        `aria-label` = "Close",
        shiny::tags$span(
          `aria-hidden` = "true"
        )
      )
    }
  )
}


## See https://getbootstrap.com/docs/5.0/components/carousel/
bs_carousel2 <- function(id, contents, interval = 4000, autostart = TRUE, wrap = TRUE, fade = FALSE) {
  autostart <- ifelse(autostart, "carousel", "false")
  wrap <- ifelse(wrap, "true", "false")
  fade <- ifelse(fade, "carousel-fade", "")
  p1 <- paste0(
    '<div id="', id, '" class="carousel slide ', fade, '" data-bs-ride="', autostart, '"',
    'data-bs-wrap="', wrap, '"><div class="carousel-inner">'
  )
  p2 <- paste0('</div>
  <button class="carousel-control-prev" type="button" data-bs-target="#', id, '" data-bs-slide="prev">
    <span class="carousel-control-prev-icon" aria-hidden="true"></span>
    <span class="visually-hidden">Previous</span>
  </button>
  <button class="carousel-control-next" type="button" data-bs-target="#', id, '" data-bs-slide="next">
    <span class="carousel-control-next-icon" aria-hidden="true"></span>
    <span class="visually-hidden">Next</span>
  </button>
  </div>')

  items <- c()
  if (length(interval) == 1) interval <- rep(interval, length(contents))
  for (i in 1:length(contents)) {
    if (i == 1) {
      div1 <- paste0('<div class="carousel-item active" data-bs-interval="', interval[i], '">')
    } else {
      div1 <- paste0('<div class="carousel-item" data-bs-interval="', interval[i], '">')
    }
    items[[i]] <- paste0(div1, contents[[i]], "</div>")
  }
  items <- paste(items, collapse = "\n")
  HTML(paste(p1, items, p2))
}


bs_carousel <- function(id, contents, img.src = NULL) {
  cbutton <- function(n) {
    shiny::tags$button(
      type = "button",
      class = ifelse(n == 1, "active", ""),
      `data-bs-target` = paste0("#", id),
      `data-bs-slide-to` = as.character(n - 1),
      `aria-current` = "true",
      `aria-label` = paste("Slide ", n)
    )
  }

  content <- ""
  citem <- function(content, img.src = NULL, active = FALSE) {
    shiny::tags$div(
      class = paste("carousel-item", ifelse(active, "active", "")),
      shiny::tags$img(
        src = img.src,
        class = "d-block w-100"
      ),
      shiny::tags$div(
        class = "carousel-caption d-none d-md-block",
        content
      )
    )
  }

  len <- length(contents)
  buttons <- lapply(1:len, cbutton)
  if (!is.null(img.src)) {
    items <- lapply(1:len, function(i) citem(contents[[i]], img.src[i], active = (i == 1)))
  } else {
    items <- lapply(1:len, function(i) citem(contents[[i]], active = (i == 1)))
  }

  shiny::div(
    id = id,
    class = "carousel slide",
    `data-bs-ride` = "carousel",
    #      shiny::div(


    #      ),
    shiny::div(
      class = "carousel-inner",
      shiny::tagList(items)
    ),
    shiny::tags$button(
      class = "carousel-control-prev",
      type = "button",
      `data-bs-target` = paste0("#", id),
      `data-bs-slide` = "prev",
      shiny::tags$span(class = "carousel-control-prev-icon", `aria-hidden` = "true"),
      shiny::tags$span(class = "visually-hidden", "Previous")
    ),
    shiny::tags$button(
      class = "carousel-control-next",
      type = "button",
      `data-bs-target` = paste0("#", id),
      `data-bs-slide` = "next",
      shiny::tags$span(class = "carousel-control-next-icon", `aria-hidden` = "true"),
      shiny::tags$span(class = "visually-hidden", "Next")
    )
  )
}

#
#

modal_carousel <- function(id, contents) {
  modalUI(
    id = id,
    title = "title",
    size = "lg",
    footer = NULL,
    bs_carousel(
      id = paste0(id, "-carousel"),
      contents = contents
    )
  )
}
