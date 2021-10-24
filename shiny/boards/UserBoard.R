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
        ),
        h3("Subscriptions"),
        shiny::div(
            id = "user-subs"
        )
    )
}

UserBoard <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE

    user <- env[["load"]][["auth"]]

    dbg("[UserBoard] >>> initializing UserBoard...")

    observeEvent(user$logged(), {
        if(!user$logged())
            return()
        
        session$sendCustomMessage(
            "get-subs",
            list(
                ns = ns(NULL)
            )
        )
    })

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
        p(
            span("Subscription level", style="color:grey;"),
            span(class = cl, tools::toTitleCase(user$level()))
        )
    })
    
    res <- list()
    return(res)
}
