
boardHeader <- function(title, info_link) {
    fillRow(
        flex=c(NA,1,NA),
        ##h2(input$nav),
        shiny::div(
            id = "navheader-current-section",
            HTML(paste0(title," &nbsp;")), 
            withTooltip( 
                shiny::actionLink(
                    inputId = info_link,
                    label="",
                    icon = shiny::icon("info-circle"),
                    style = "color: #ccc;"
                ),
                "Show information and tutorial about this board"
            )
        ),        
        shiny::br(),
        ##shiny::div(shiny::textOutput("current_dataset"), class='current-dataset')
        shiny::div("Data set name selected", id='navheader-current-dataset')        
    )
}
