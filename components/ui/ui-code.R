##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## textInput <- function(inputId, label, value = "") {
myTextInput <- function(inputId, label, value = "") {
  #
  shiny::tagList(
    shiny::tags$label(label, `for` = inputId),
    shiny::tags$input(
      id = inputId, type = "text", value = value,
      class = "myTextInput form-control shiny-bound-input"
    )
  )
}

## shiny::shinyUI(
##     shiny::basicPage(
##         code
#
#
#
##     ))

## from https://gist.github.com/xiaodaigh/7150112
code.textInput <- shiny::HTML(" <script> var myTextInputBinding = new Shiny.InputBinding();
            $.extend(myTextInputBinding, {
            find: function(scope) {
            return $(scope).find('.myTextInput');
            },
            getId: function(el) {
            //return InputBinding.prototype.getId.call(this, el) || el.name;
            return $(el).attr('id')
            },
            getValue: function(el) {
            return el.value;
            },
            setValue: function(el, value) {
            el.value = value;
            },
            subscribe: function(el, callback) {
             $(el).on('keyup.textInputBinding input.textInputBinding', function(event) {
               if(event.keyCode == 13) { //if enter
                callback()
               }
             });
            $(el).on('focusout.myTextInputBinding', function(event) { // on losing focus
              callback();
            });
            },
            unsubscribe: function(el) {
            $(el).off('.myTextInputBinding');
            },
            receiveMessage: function(el, data) {
            if (data.hasOwnProperty('value'))
            this.setValue(el, data.value);

            if (data.hasOwnProperty('label'))
            $(el).parent().find('label[for=' + el.id + ']').text(data.label);

            $(el).trigger('change');
            },
            getState: function(el) {
            return {
            label: $(el).parent().find('label[for=' + el.id + ']').text(),
            value: el.value
            };
            },
            getRatePolicy: function() {
            return {
            policy: 'debounce',
            delay: 250
            };
            }
            });
            Shiny.inputBindings.register(myTextInputBinding, 'shiny.myTextInput');</script>")

inlineTextInput <- function(inputId, label, value = "") {
  shiny::div(
    style = "display:inline-block",
    shiny::tags$label(label, "for" = inputId),
    shiny::tags$input(id = inputId, type = "text", value = value, class = "input-small")
  )
}

inlineSelectInput <- function(inputId, label, value = "") {
  shiny::div(
    style = "display:inline-block",
    shiny::tags$label(label, "for" = inputId),
    shiny::tags$input(id = inputId, type = "select", value = value)
  )
}




# Copied code from shiny::fileInput but changed the class of the button
# it creates to btn-outline-primary.
# Maybe on future versions of that function they will add
# a `class` argument so we can submit our own class to that button.
fileInput2 <- function(
    inputId, label, multiple = FALSE, accept = NULL, width = NULL,
    buttonLabel = "Browse...", placeholder = "No file selected",
    capture = NULL, buttonClass = "btn-outline-primary") {
  restoredValue <- restoreInput(id = inputId, default = NULL)
  if (!is.null(restoredValue) && !is.data.frame(restoredValue)) {
    warning("Restored value for ", inputId, " has incorrect format.")
    restoredValue <- NULL
  }
  if (!is.null(restoredValue)) {
    restoredValue <- toJSON(restoredValue, strict_atomic = FALSE)
  }
  inputTag <- tags$input(
    id = inputId, name = inputId, type = "file",
    style = "position: absolute !important; top: -99999px !important; left: -99999px !important;",
    `data-restore` = restoredValue
  )
  if (multiple) {
    inputTag$attribs$multiple <- "multiple"
  }
  if (length(accept) > 0) {
    inputTag$attribs$accept <- paste(accept, collapse = ",")
  }
  if (!is.null(capture)) {
    inputTag$attribs$capture <- capture
  }
  shiny::tags$div(
    class = "form-group shiny-input-container",
    style = htmltools::css(width = htmltools::validateCssUnit(width)),
    shiny:::shinyInputLabel(inputId, label), div(
      class = "input-group",
      shiny::tags$label(
        class = "input-group-btn input-group-prepend",
        span(
          class = paste("btn btn-default", buttonClass), buttonLabel,
          inputTag
        )
      ), tags$input(
        type = "text", class = "form-control",
        placeholder = placeholder, readonly = "readonly"
      )
    ),
    shiny::tags$div(
      id = paste(inputId, "_progress", sep = ""),
      class = "progress active shiny-file-input-progress",
      shiny::tags$div(class = "progress-bar")
    )
  )
}

selector_switch <- function(
    class = NULL,
    label = "Text to appear in Switch",
    is.checked) {
  tags$div(
    class = "form-check form-switch",
    tags$input(
      class = paste("form-check-input", class),
      type = "checkbox",
      id = ifelse(is.checked, "flexSwitchCheckChecked", "flexSwitchCheckDefault"),
      role = "switch"
    ),
    tags$label(
      class = "form-check-label",
      `for` = "flexSwitchCheckDefault",
      label
    )
  )
}

loading_spinner <- function(text = "Loading...") {
  shiny::tags$div(
    id = "spinner-container", # Add an ID to the spinner container
    class = "spinner-container",
    shiny::tags$div( # Wrap the spinner and text in an additional div
      class = "spinner-wrapper",
      shiny::tags$div(
        class = "spinner-border text-primary", role = "status",
        style = "font-size:0.4em;height:4em;width:4em;",
        shiny::tags$span(class = "visually-hidden", "Loading...")
      ),
      shiny::tags$p(class = "spinner-text", text) # Add a class to the text element
    )
  )
}
