
##textInput <- function(inputId, label, value = "") {
myTextInput <- function(inputId, label, value = "") {
    ##singleton(shiny::tags$head(shiny::tags$script(src = "/temp/mytextinput.js"))),
    shiny::tagList(shiny::tags$label(label, `for` = inputId),
            shiny::tags$input(id = inputId, type = "text", value = value,
                       class="myTextInput form-control shiny-bound-input"))
}

## shiny::shinyUI(
##     shiny::basicPage(
##         code
##        ,myTextInput("myTextInput","My text input","On enter or focus out")
##        ,shiny::textOutput("meh")
##        ,shiny::HTML('<script src="https://gist.github.com/xiaodaigh/7150112.js"></script>')
##     ))

## from https://gist.github.com/xiaodaigh/7150112
code.textInput = shiny::HTML(" <script> var myTextInputBinding = new Shiny.InputBinding();
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

inlineTextInput <- function (inputId, label, value = "")
{
    shiny::div(style="display:inline-block",
        shiny::tags$label(label, "for"=inputId),
        shiny::tags$input(id=inputId, type="text", value=value, class="input-small"))
}

inlineSelectInput <- function (inputId, label, value = "")
{
    shiny::div(style="display:inline-block",
        shiny::tags$label(label, "for"=inputId),
        shiny::tags$input(id=inputId, type="select", value=value))
}


# This JS code is responsible for activating and desactivating the class `active`,
# defined on `scss/components/_dropdownmenu.scss` when clicking a dropdown button.
# It also has logic implemented so that when clicking a different dropdown item,
# the active one gets desactivated.
jsCode <- function(id){
  paste0(
    "$('#", id, "').on('click', function(){
    if($('#", id, "').hasClass('active')){
         $('#", id,"').toggleClass('active');
    return 0;
      }
      $('.dropdown-button.active').toggleClass('active');
      if(!$('#", id, "').hasClass('active')){
         $('#", id,"').toggleClass('active');
      }

})"
  )
}

DropdowMenu <- function (..., size = "default", status = "default", icon = NULL, margin = "10px")
{
  id <- bigdash:::make_id()
  tags$div(tags$a(class = paste0(
    "btn btn-", status, " action-button dropdown-button btn-circle", if(size == "default"){" "}else{paste0("-", size, " ")},
    if (size != "default") paste0("btn-", size)
  ),
  id = id, type = "button", `data-bs-toggle` = "dropdown",
  `aria-expanded` = "false", list(icon)),
  tags$div(class = "dropdown-menu",
           `aria-labelledby` = id,
           tryCatch({
             lapply(list(...), function(x){
               htmltools::tags$li(
                 x,
                 style = paste0("margin-left: ",
                                htmltools::validateCssUnit(margin), "; margin-right: ",
                                htmltools::validateCssUnit(margin), ";")
               )
             })
           }, error = function(w){})
  ),
  tags$script(HTML(jsCode(id)))
  )
}