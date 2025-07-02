fileInputArea <- function(inputId, label, multiple = FALSE, accept = NULL,
                          width = NULL,
                          buttonLabel = "Drag your file here or click to browse",
                          placeholder = "No file selected",
                          fileBrowser = FALSE) {
  restoredValue <- restoreInput(id = inputId, default = NULL)

  # Catch potential edge case - ensure that it's either NULL or a data frame.
  if (!is.null(restoredValue) && !is.data.frame(restoredValue)) {
    warning("Restored value for ", inputId, " has incorrect format.")
    restoredValue <- NULL
  }

  if (!is.null(restoredValue)) {
    restoredValue <- toJSON(restoredValue, strict_atomic = FALSE)
  }

  inputTag <- tags$input(
    id = inputId,
    name = inputId,
    type = "file",
    # Don't use "display: none;" style, which causes keyboard accessibility issue; instead use the following workaround: https://css-tricks.com/places-its-tempting-to-use-display-none-but-dont/
    style = "position: absolute !important; top: -99999px !important; left: -99999px !important;",
    `data-restore` = restoredValue
  )

  if (multiple) {
    inputTag$attribs$multiple <- "multiple"
  }
  if (length(accept) > 0) {
    inputTag$attribs$accept <- paste(accept, collapse = ",")
  }

  if (fileBrowser) {
    return(shinyfilebrowser::file_browser_ui(inputId))
  }

  div(
    class = "form-group", # shiny-input-container w-100",
    style = htmltools::css(width = htmltools::validateCssUnit(width), margin = "auto"),
    shiny:::shinyInputLabel(inputId, ""),
    div(
      class = "input-group mb-3",
      # input-group-prepend is for bootstrap 4 compat
      tags$label(
        class = "input-group-btn input-group-prepend w-100",
        span(
          class = "btn btn-area w-100", inputTag,
          ##          div(tags$image(src = fileInputArea_icon_encoded, width = "80px;"), style = "margin-top: 2rem;"),
          div(icon("upload"), style = "font-size: 80px;"),
          div(p(label), style = "font-size: 1.2rem; font-weight: 700; padding-top: 2rem;"),
          div(p(buttonLabel), style = "font-size: 1rem; font-weight: 400; margin: 0.7rem;"),
          div("browse", style = "font-size:0.92rem; color:white; padding:6px 6px; background-color:#3181DE; width:75px; margin:auto;border-radius:4px;") ## fake button
        )
      )
    ),
    tags$div(
      id = paste0(inputId, "_progress"),
      class = "progress active shiny-file-input-progress",
      tags$div(class = "progress-bar")
    )
  )
}

## Icon from <https://icons.getbootstrap.com/icons/upload/>
fileInputArea_icon_file <- tempfile(fileext = ".svg")
writeLines('
<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="#495057" class="bi bi-upload" viewBox="0 0 16 16">
  <path d="M.5 9.9a.5.5 0 0 1 .5.5v2.5a1 1 0 0 0 1 1h12a1 1 0 0 0 1-1v-2.5a.5.5 0 0 1 1 0v2.5a2 2 0 0 1-2 2H2a2 2 0 0 1-2-2v-2.5a.5.5 0 0 1 .5-.5z"/>
  <path d="M7.646 1.146a.5.5 0 0 1 .708 0l3 3a.5.5 0 0 1-.708.708L8.5 2.707V11.5a.5.5 0 0 1-1 0V2.707L5.354 4.854a.5.5 0 1 1-.708-.708l3-3z"/>
</svg>',
  con = fileInputArea_icon_file
)
fileInputArea_icon_encoded <- xfun::base64_uri(fileInputArea_icon_file) ## can we not pass a string??
