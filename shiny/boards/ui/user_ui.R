##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

UserInputs <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
               HTML("<h3>User Settings</h3><br><br>"),
               shiny::uiOutput(ns("description"))
               ## shiny::uiOutput(ns("inputsUI"))
           )
}

UserUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
               height = 750,
               shiny::tabsetPanel(
                          id = ns("tabs"),
                          shiny::tabPanel("User settings",uiOutput(ns("userinfo_UI")))
                          ## shiny::tabPanel("Visitors map",uiOutput(ns("usersmap_UI")))
                          ## shiny::tabPanel("Community forum",uiOutput(ns("forum_UI")))
                      )
           )
}