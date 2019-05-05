
##textInput <- function(inputId, label, value = "") {
myTextInput <- function(inputId, label, value = "") {
    ##singleton(tags$head(tags$script(src = "/temp/mytextinput.js"))),
    tagList(tags$label(label, `for` = inputId),
            tags$input(id = inputId, type = "text", value = value,
                       class="myTextInput form-control shiny-bound-input"))
}

## shinyUI(
##     basicPage(
##         code
##        ,myTextInput("myTextInput","My text input","On enter or focus out")
##        ,textOutput("meh")
##        ,HTML('<script src="https://gist.github.com/xiaodaigh/7150112.js"></script>')
##     ))

## from https://gist.github.com/xiaodaigh/7150112
code.textInput = HTML(" <script> var myTextInputBinding = new Shiny.InputBinding();
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

textInputRow <- function (inputId, label, value = "")
{
    div(style="display:inline-block",
        tags$label(label, "for"=inputId),
        tags$input(id=inputId, type="text", value=value, class="input-small"))
}

selectInputRow <- function (inputId, label, value = "")
{
    div(style="display:inline-block",
        tags$label(label, "for"=inputId),
        tags$input(id=inputId, type="select", value=value))
}

if(0) {
    runApp(list(
        ui = bootstrapPage(
            textInputRow(inputId="xlimitsmin", label="x-min", value = 0.0),
            textInputRow(inputId="xlimitsmax", label="x-max", value = 0.5)
        ),
        server = function(input, output) {}
    ))
}
