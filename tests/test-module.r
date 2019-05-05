library(shiny)
source("module-linkedScatter.R")

ui <- fixedPage(
    h2("Module example"),
    selectInput("x1","x1",c("cty","drv","cyl","manufacturer")),
    selectInput("x2","x2",c("cty","drv","cyl","manufacturer")),
    linkedScatterUI("scatters"),
    textOutput("summary")
)

server <- function(input, output, session) {
    ##  df <- callModule(linkedScatter, "scatters", reactive(mpg),
    ##    left = reactive(c("cty", "hwy")),
    ##    right = reactive(c("drv", "hwy"))
    ##  )
    df <- callModule(linkedScatter, "scatters", mpg,
                     ##left = c("cty", "hwy"), right = c("drv", "hwy")
                     left  = reactive(c(input$x1, "hwy")),
                     right = reactive(c(input$x2, "hwy"))
                     )

    output$summary <- renderText({
        sprintf("%d observation(s) selected", nrow(dplyr::filter(df(), selected_)))
    })
}

shinyApp(ui, server)
