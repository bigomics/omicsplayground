##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##


QUESTIONS.ANSWERS <- list(
    c("Would you recommend Omics Playground to your colleagues?","no|maybe|yes|sure|definitely!"),
    c("What is your favorite fruit?","banana|apple|pear"),
    c("Would you like to receive our newsletter?","yes|no"),
    c("Your email please?","<text>"),
    c("Are you a student?","yes|no")
)

QuestionBoard_UI <- function(id) {
    ##ns <- NS(id)  ## namespace
}

QuestionBoard <- function(input, output, session, lapse=5)
{
    ns <- session$ns ## NAMESPACE
    
    ##qa = QUESTIONS.ANSWERS[[1]]
    cur_question <- reactiveVal("")
    
    showQuestion <- function(qa=NULL) {
        dbg("[QuestionBoard:showQuestion] reacted")
        if(is.null(qa)) qa <- sample(QUESTIONS.ANSWERS,1)[[1]]
        question <- qa[1]
        cur_question(qa[1])
        choices <- strsplit(qa[2],split="\\|")[[1]]
        if(choices[1]=="<text>") {
            answer.tags <- textInput(ns("answer"),NULL, value="")
        } else {
            answer.tags <- radioButtons(ns("answer"),NULL, choices=choices)
        }

        showModal( modalDialog(
            title = question,
            answer.tags,
            footer = tagList(
                ##modalButton("Cancel"),
                actionButton(ns("question_submit"),"Continue", icon=NULL)
            ),
            size = "m"
        ))
    }

    nmesg = 0
    observe({
        if(lapse>0 && nmesg >= 1) {
            showQuestion(qa=NULL)
        }
        if(lapse>0) invalidateLater(1000*60*lapse)  ## every 10 minutes??
        nmesg <<- nmesg + 1
    })

    observeEvent( input$question_submit, {
        dbg("[QuestionBoard$:observeEvent] input$question_submit")
        question <- cur_question()
        answer <- input$answer
        qa <- data.frame(time=date(), question="question", answer="answer")
        qa <- data.frame(time=date(), question=question, answer=answer)
        dbg("[QuestionBoard$:observeEvent] updating answer.csv file...")
        write.table( qa, file="answers.csv", sep=",", append=TRUE,
                    row.names=FALSE, col.names=FALSE)
        removeModal()        
    }) 
    
    
}

