##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


header_infotext <- "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."

boardHeader <- function(title, info_link) {
  div(
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
                    ## icon = shiny::icon("info-circle"),
                    icon = shiny::icon("youtube"),                    
                    style = "color: #ccc;"
                ),
                "Video tutorial about this board"
            )
        ),
        div(class = "card-footer-switch", 
            selector_switch(
                class = 'card-footer-checked', 
                label = "show captions",
                is.checked = FALSE
                )
            )
    )
##    bs_alert(shiny::HTML(header_infotext))
  )
}
