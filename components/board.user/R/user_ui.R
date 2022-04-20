##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

UserInputs <- function(id) {
    ns <- shiny::NS(id)
    bigdash::tabSettings(
        shiny::uiOutput(ns("description"))
    )
}

UserUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
               height = 750,
               shiny::tabsetPanel(
                          id = ns("tabs"),
                          shiny::tabPanel("User settings",
                          fillRow(
                            flex=c(0.8,0.2,1,0.2,1),
                            tagList(
                                shiny::h4("News"),            
                                shiny::htmlOutput(ns("news"))
                            ),br(),
                            tagList(
                                shiny::h4("Personal"),
                                uiOutput(ns("plan")),                    
                                shiny::tableOutput(ns("userdata"))
                            ),br(),
                            tagList(
                                shiny::h4("Settings"),            
                                shinyWidgets::prettySwitch(ns("enable_beta"),"enable beta features")
                            )
                        ))
                        # Currently not used Stefan 22.03.22
                        # shiny::tabPanel("Visitors map",
                        #     shiny::fillCol(
                        #         height = 600,
                        #         shiny::fillRow(
                        #             flex = c(1,4.5),
                        #             shiny::wellPanel( shiny::uiOutput(ns("usersmapInfo"))),
                        #             plotWidget(ns("usersmap"))
                        #         )
                        #     )
                        # )
                        #shiny::tabPanel("Community forum",uiOutput(ns("forum_UI")))
                      )
           )
}