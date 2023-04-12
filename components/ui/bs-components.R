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


bs_carousel <- function(id, contents, img.src=NULL) {
    
  cbutton <- function(n) {
    shiny::tags$button(
      type = "button",
      class = ifelse(n==1, "active",""),
      `data-bs-target` = paste0("#",id),
      `data-bs-slide-to` = as.character(n-1),
      `aria-current` = "true",
      `aria-label` = paste("Slide ",n)
      )
  }

  content=""
  citem <- function(content, img.src=NULL, active=FALSE) {
    shiny::tags$div(
      class = paste("carousel-item",ifelse(active,"active","")),
      shiny::tags$img(
         src = img.src,
         class="d-block w-100"
      ),
      shiny::tags$div(
        class = "carousel-caption d-none d-md-block",
        content
      )
    )
  }

  len <- length(contents)
  buttons <- lapply(1:len, cbutton)
  if(!is.null(img.src)) {
      items   <- lapply(1:len, function(i) citem(contents[[i]], img.src[i], active=(i==1)))
  } else {
      items   <- lapply(1:len, function(i) citem(contents[[i]], active=(i==1)))
  }     
  
  shiny::div(
      id = id,
      class = "carousel slide",
      `data-bs-ride` = "carousel",
      shiny::div(
          class="carousel-indicators",
          shiny::tagList(buttons)
      ),
      shiny::div(
          class="carousel-inner",
          shiny::tagList(items)
      ),
      shiny::tags$button(
        class = "carousel-control-prev",
        type = "button",
        `data-bs-target` = paste0("#",id),
        `data-bs-slide` = "prev",
        shiny::tags$span( class = "carousel-control-prev-icon", `aria-hidden` = "true"),
        shiny::tags$span( class = "visually-hidden", "Previous")
      ),
      shiny::tags$button(
        class = "carousel-control-next",
        type = "button",
        `data-bs-target` = paste0("#",id),
        `data-bs-slide` = "next",
        shiny::tags$span( class = "carousel-control-next-icon", `aria-hidden` = "true"),
        shiny::tags$span( class = "visually-hidden", "Next" )
      )
  )
}

##contents=list("hello","world");id="id"
##bs_carousel(id=id, contents) 

modal_carousel <- function(id, contents) {
    modalUI(
        id = id,
        title = "title",
        size = "lg",
        footer = NULL,
        bs_carousel(
            id = paste0(id,"-carousel"),
            contents = contents
        )
    )
}
