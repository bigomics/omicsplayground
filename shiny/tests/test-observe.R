library(shiny)

testModuleUI <- function(id) {
    ns <- NS(id)
    tagList(
        actionButton(ns("button2"), "button2"),
        verbatimTextOutput(ns("out"))
    )
}

testModule <- function(input, output, session) {
    count2 <- eventReactive(input$button2, {
        input$button2
    })    
    output$out <- renderText(count2())
}

ui = fillPage(
    actionButton("button1",label="button1"),
    verbatimTextOutput("out"),
    testModuleUI("test1")
)

server <- function(input, output, session) {
    callModule(testModule, "test1")

    count1 <- eventReactive(input$button1, {
        input$button1
    })    
    output$out <- renderText(count1())
}

shinyApp(ui, server)
