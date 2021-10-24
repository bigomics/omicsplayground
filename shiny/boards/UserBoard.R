##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing UserBoard")

UserInputs <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::uiOutput(ns("description"))
    )
}

UserUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace

    shiny::div(
        style = "padding-left:1rem!important;",
        shiny::h1("User Profile"),
        shiny::div(
            class = "row",
            shiny::div(
                class = "col-md-3",
                uiOutput(ns("username")),
                uiOutput(ns("plan"))
            )
        )
    )
}

UserBoard <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE

    user <- env[["load"]][["auth"]]

    dbg("[UserBoard] >>> initializing UserBoard...")

    output$username <- renderUI({
        div(
            h3(user$name()),
            p(user$email())
        )
    })

    output$plan <- renderUI({
        plan_class <- "info"
        if(user$level() == "premium")
            plan_class <- "success"

        cl <- sprintf("badge badge-%s", plan_class)
        span(class = cl, tools::toTitleCase(user$level()))
    })
    
    output$userdata <- renderTable({
        dbg("[UserBoard::userdata]  renderDataTable")
        values <- c(
            name   = user$name(),
            email  = user$email(),
            plan   = user$level(),
            ##logged = user$logged(),
            limit  = paste(user$limit(),collapse=';')
        )
        values[which(values=="")] <- "(not set)"
        df <- data.frame(field=names(values), value=values)
    })
        
    res <- list()
    return(res)
}
