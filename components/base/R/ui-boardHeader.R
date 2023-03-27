##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

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
        div(class = "card-selector-header", 
            selector_default(
                class = 'card-footer-checked', 
                label = "Show captions"))
    )
}
