library(shiny)

counterButton <- function(id, label = "Counter") {
    ns <- NS(id)
    tagList(
        actionButton(ns("button"), label = label),
        verbatimTextOutput(ns("out"))
    )
}

counter <- function(input, output, session, start=0) {
    count <- reactiveVal(start)
    observeEvent(input$button, {
        count(count() + 1)
    })
    output$out <- renderText({
        count()
    })
    count
}

ui <- fluidPage(
    counterButton("counter1", "Counter #1"),
    counterButton("counter2", "Counter #2")
)

server <- function(input, output, session) {
    callModule(counter, "counter1", start=0)
    callModule(counter, "counter2", start=10)
}

shinyApp(ui, server)
