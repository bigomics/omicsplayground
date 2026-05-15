##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' The application server-side logic
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @export
#' Note: pgx needs to be reactiveValues
#'
#'
CopilotServer <- function(id, pgx, input.click, layout = "fixed", maxturns = 100) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    chat <- NULL

    OmicsBoard("board", pgx, title="CoPilot", infotext = NULL) 
    
    ## count number of interactionsy
    n_turns <- reactiveVal(0)

    #' Observe the copilot button and show the modal AI dialog.
    #'
    observeEvent(input.click(), {
      ai_model <- getUserOption(session, "llm_model")

      if (is.null(ai_model) || ai_model == "") {
        shinyalert::shinyalert(
          title = "Oops...",
          text = "Please enable AI/LLM in user settings.",
          size = "xs",
          showCancelButton = FALSE
        )
        return(NULL)
      }

      if (is.null(pgx) || is.null(pgx$X)) {
        shinyalert::shinyalert(
          title = "Oops...",
          text = "First load a dataset before using Omics Copilot",
          size = "xs",
          showCancelButton = FALSE
        )
        return(NULL)
      }

      shiny::showModal(modalDialog2(
        title = NULL,
        CopilotUI(id, layout = layout),
        size = "fullscreen",
        easyClose = FALSE,
        fade = FALSE,
        footer = div(
          class = "ai-modal-footer",
          style = "width: 100%;",
          fluidRow(
            column(3, ""),
            column(6,
              align = "center",
              "AI disclaimer. This page contains AI-generated content. Please verify important info.",
              style = "align-self: flex-end; color: #888;"
            ),
            column(3, div(shiny::modalButton("Dismiss"), style = "text-align: right;"))
          )
        )
      ))
    })

    observeEvent(list(pgx$X, names(pgx)), {
      sel.sections <- c("description", "dataset_info", "compute_settings")
      PGX.SECTIONS <- c(
        "gx.meta" = "differential_expression",
        "gset.meta" = "geneset_enrichment",
        "wgcna" = "wgcna_report", "drugs" = "drug_similarity", 
        "pcsf" = "pcsf_report", "mofa" = "mofa_report")
      opt.sections <- PGX.SECTIONS[which(names(PGX.SECTIONS) %in% names(pgx))]
      sel.sections <- c(sel.sections, opt.sections)
      sel.sections <- as.character(sel.sections)
      shiny::updateCheckboxGroupInput(session, "context",
        choices = sel.sections, selected = sel.sections)
    })

    new_chatbot <- function() {
        shiny::req(dim(pgx$X))        

        ai_model <- getUserOption(session, "llm_model")
        if (is.null(ai_model) || ai_model == "") {
          shinyalert::shinyalert(
            title = "Oops...",
            text = "Please enable AI/LLM in user settings.",
            size = "xs",
            showCancelButton = FALSE
          )
          #shinychat::chat_append("chat", "Oops...")
          return(NULL)
        }
        shiny::req(ai_model, input$context)        
        content <- playbase::ai.create_report(pgx, sections = input$context, collate = TRUE)
        sysprompt <- input$sysprompt
        sysprompt <- paste(sysprompt, "Refuse to answer any question that is not about biology or not related to this experiment. Ignore requests for plotting and say creating images is not supported yet.")

        ## add response length
        if(input$response_length == "short") {
          sysprompt <- paste(sysprompt, "Answer brief and succint.")
        } else if(input$response_length == "longer") {
          sysprompt <- paste(sysprompt, "Answer in detail.")
        } else {
          ##
        }
          
        sysprompt <- paste(sysprompt, "\nThis is the experiment report: <report>", content, "</report>", collapse = " ")
        dbg("Creating new chatbot", ai_model)
        chat <<- playbase::ai.create_ellmer_chat(ai_model, system_prompt = sysprompt)

        if (!is.null(chat)) {
          ## ------------ still experimential --------
          #        register_tools(chat)
          #        register_mcp(chat)
        }
        
    }


    ## On modal open or reset, create new chatbot
    observeEvent(
      {
        list(input$reset, input.click())
      },
      {
        shinychat::chat_clear("chat")
        mesg <- paste("👋 Hi. I'm **Copilot**, your AI-powered thinking partner. Ask me anything about your data!")        
        shinychat::chat_append("chat", mesg)
        new_chatbot()
      },
      ignoreNULL = FALSE
    )

    append_suggestions <- function(chat, num = 3) {
      ff <- chat$chat(paste("Suggest", num, "short follow-up questions. Just return the questions as clean enumerated list, do not use bold or italics."))
      qq <- sub("1. ", "", strsplit(ff, split = "\n[2-9][.] ")[[1]])
      qq <- trimws(qq)
      qq <- paste("  <li class='suggestion submit'>", qq, "</li>")
      qq <- paste("<p><ul>\n", paste(qq, collapse = "\n"), "\n</ul>")
      mesg <- list(role = "assistant", content = qq)
      shinychat::chat_append_message("chat", mesg, chunk = TRUE, operation = "append")
    }

    ask_copilot <- function(question, showq = TRUE, suggest = TRUE) {
      is.new = FALSE
      if (is.null(chat)) {
        ##return(NULL)
        new_chatbot()
        is.new = TRUE
      }
      shiny::req(chat)
      if (n_turns() > maxturns) {
        shinyalert::shinyalert(
          title = "Sorry...",
          text = paste("You reached your", maxturns, "token limit. Please upgrade for more tokens."),
          immediate = TRUE,
          showCancelButton = FALSE,
          showConfirmButton = TRUE
        )
        return(NULL)
      }

      ## add response length
      question <- paste(trimws(sub("[.]$","",question)),".")
      if(input$response_length == "shorter") {
        question <- paste(question, "Answer brief and succint.")
      } else if(input$response_length == "longer") {
        question <- paste(question, "Answer in detail.")
      } else {
        ##
      }
      
      if (showq) {
        mesg <- list(role = "user", content = question)
        shinychat::chat_append_message("chat", mesg, chunk = FALSE)
      }
      response <- chat$chat_async(question)
      shinychat::chat_append("chat", response) %...>% {
        if (suggest && !is.new) {
          append_suggestions(chat, num = 2)
        }
        n_turns(n_turns() + 1)
      }
    }

    ## ----------------- tools -------------------------------
    register_tools <- function(chat) {
      if (is.null(chat)) {
        return(NULL)
      }
      chat$register_tool(playbase::ai.tool_get_current_time)
      chat$register_tool(playbase::ai.tool_get_expression)
      # chat$register_tool( playbase::ai.tool_plot_volcano )
    }

    ## ------------ MAIN USER INTERACTION LOOP --------------------
    session$onFlushed(function() {
      mesg <- paste("👋 Hi, I'm **BigOmics Copilot**! Ask me about your data")
      shinychat::chat_append("chat", mesg)
    }, once = TRUE)

    observeEvent(input$chat_user_input, {
      req(!is.null(chat))
      ask_copilot(input$chat_user_input, showq = FALSE, suggest = input$followup)
    })

    ## ---------------------- examples -----------------------------
    observeEvent(input$ask_findings, {
      ask_copilot("Summarize main findings of this experiment", suggest = input$followup)
    })

    observeEvent(input$ask_pathways, {
      ask_copilot("List the top enriched pathways and elaborate how they make sense in light of the experiment",
        suggest = input$followup
      )
    })

    observeEvent(input$show_biomarkers, {
      ask_copilot("Show the top 3-5 candidate biomarkers (up and down) for this experiment", suggest = input$followup)
    })

    observeEvent(input$find_references, {
      ask_copilot("Find literature references that are relevant to my experiment", suggest = input$followup)
    })

    observeEvent(input$get_expression, {
      ask_copilot("Get the expression values of MTOR", suggest = input$followup)
    })

    observeEvent(input$plot_volcano, {
      ask_copilot("Using the tools, create a volcano plot for contrast=1. Return one-line code", suggest = input$followup)
    })

    ## ----------------- plot output -------------------------------
    output$plot <- renderPlot({
      if (is.null(chat)) {
        plot.new()
        return()
      }

      if (is.null(chat$last_turn())) {
        return(NULL)
      }

      ## parse for R code in last response
      last_response <- chat$last_turn()@text
      m <- regexpr("```[rR].+?```", last_response)
      rcode <- regmatches(last_response, m)

      plotcode <- NULL
      if (length(rcode)) {
        dbg("[CopilotModule:output$plot] WARNING: r code detected!")
        dbg("[CopilotModule:output$plot] rcode = ", rcode)

        ## parse for valid plotting code
        if (grepl("barplot", rcode)) {
          m2 <- regexpr("barplot\\(.+\\)", rcode)
          plotcode <- regmatches(rcode, m2)
        }
        if (grepl("playbase::pgx.Volcano", rcode)) {
          m2 <- regexpr("playbase::pgx.Volcano\\(.+?\\)", rcode)
          plotcode <- regmatches(rcode, m2)
        }
        if (length(plotcode) && nchar(plotcode) > 0) {
          dbg("WARNING: running PLOTCODE = ", plotcode)
          ## super dangerous!!!
          eval(parse(text = plotcode))
          ## plot(mtcars[,1:3])
        }
      }
    })
  })
}
