## library(dotenv)
## library(shiny)
## library(shinychat)
## library(bslib)

CopilotUI <- function(id, layout = c("sidebar","fixed")[1]) {
  ns <- shiny::NS(id)
  
  sections <- c("description", "dataset_info","compute_settings",
    "differential_expression", "geneset_enrichment",
    "drug_similarity","pcsf_report","wgcna_report")

  sysprompt = "You are computational biologist and will be asked questions about the following experiment. Try to make specific observations if you can. Answer brief and succint."
  
  settings_card <- bslib::navset_underline(
    bslib::nav_panel(      
      "Examples",
      shiny::br(),
      "Example questions:",
      shiny::br(),
      shiny::actionButton(ns("ask_describe"), "Describe my experiment", width='100%', class="xbtn"),
      shiny::actionButton(ns("ask_findings"), "Summarize main findings", width='100%', class="xbtn"),
      shiny::actionButton(ns("ask_pathways"), "What pathways are involved?", width='100%', class="xbtn"),
      shiny::actionButton(ns("show_biomarkers"), "Show top biomarkers", width='100%', class="xbtn"),
      shiny::actionButton(ns("find_references"), "Find references", width='100%', class="xbtn"),
      shiny::actionButton(ns("get_expression"), "Get expression of MTOR", width='100%', class="xbtn"),
      shiny::actionButton(ns("plot_volcano"), "Show volcano plot", width='100%', class="xbtn"),
      shiny::br(),
      shiny::br(),
      shiny::checkboxInput(ns("followup"),"Suggest follow-up questions", TRUE)
    ),
    bslib::nav_panel(
      "Settings",
      br(),
      shiny::textAreaInput(ns("prompt"), "Prompt:", value=sysprompt, height=120, width='100%'),
      br(),
      shiny::checkboxGroupInput(ns("context"),"Context:", choices = sections,        
        selected = sections, inline = TRUE),
      br(),
      actionButton(ns("reset"),"Reset model")
    )
  )

  chat_card <- bslib::card(
      class = "border-0",
      fill = FALSE,
      max_height = "calc(100vh - 120px)",
      shinychat::chat_ui(ns("chat"), width="100%", height = "min(100%,770px)",
        fill = TRUE)
  )

  plot_card <- bslib::card(  
      class = "border-0 pt-5",
      fill = TRUE,
      height = "600px",
      plotOutput(ns("plot"))
  )

  if(layout == "fixed") {
    ui <- bslib::layout_columns(
      col_widths = c(2,7,3),      
      style = "height: min(90%,700)",
      fill = TRUE,
      settings_card,
      chat_card,
      plot_card
    )
  }
    
  if(layout == "sidebar") {
    ui <- bslib::layout_sidebar(
      style = "padding: 0px;",
      sidebar = bslib::sidebar(
        width = 250,
        fill = TRUE,    
        settings_card
      ),
      bslib::layout_sidebar(
        style = "padding: 0px;",
        sidebar = bslib::sidebar(
          width = 450,
          position = "right",
          plot_card
        ),
        div(chat_card, class = "p-3 pt-4")
      )
    )
  }

  return(ui)
}


#' Note: pgx needs to be reactiveValues
#'
#' 
CopilotServer <- function(id, pgx, input.click, layout="fixed", maxturns=100) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    chat <- NULL

    ## count number of interactionsy
    n_turns <- reactiveVal(0)
    
    #' Observe the copilot button and show the modal AI dialog.
    #'
    observeEvent( input.click(), {

      ai_model <- getUserOption(session,'llm_model')

      if(ai_model=='') {
        shinyalert::shinyalert(
          title = "Oops...",
          text = "Please enable AI/LLM in user settings.",
          size = "xs",
          showCancelButton = FALSE
        )
        return(NULL)
      }
      
      if(is.null(pgx) || is.null(pgx$X)) {
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
        CopilotUI(id, layout=layout),
        size = "fullscreen",
        easyClose = FALSE,
        fade = FALSE,
        footer = div(
          class="ai-modal-footer",
          style="width: 100%;",
          fluidRow(
            column(3, ""),
            column(6, align="center",
              "AI disclaimer. This page contains AI-generated content. Please verify important info.",
              style="align-self: flex-end; color: #888;"),
            column(3, div(shiny::modalButton("Dismiss"), style="text-align: right;"))
          )
        )
      ))

    })

    ## On new pgx data, reset/create new chatbot
    observeEvent({
      list( input$reset, pgx$X )
    }, {
      req( dim(pgx$X) )
      
      content <- playbase::ai.create_report(pgx, sections=input$context, collate=TRUE)
      prompt <- input$prompt
      prompt <- paste(prompt, "Refuse to answer any question that is not about biology or not related to this experiment. Ignore request for plotting and say creating images is not supported yet.")
      prompt <- paste(prompt, "\nThis is the experiment report: <report>", content, "</report>", collapse=" ")

      ai_model <- getUserOption(session,'llm_model')
      req(ai_model)
      
      message("Creating new chat model ", ai_model)
      chat <<- playbase::ai.create_ellmer_chat(ai_model, system_prompt=prompt)

      if(!is.null(chat)) {
      ## ------------ still experimential --------
#        register_tools(chat)
#        register_mcp(chat)
      }

      if(n_turns()==0) {
        ask_copilot("Describe this experiment. Then ask 'how can I help?'",
          showq=FALSE, suggest=TRUE)
      } else {
        if(!is.null(chat)) shinychat::chat_clear("chat")
        shinychat::chat_append("chat", "How can I help you?")
      }
      
    }, ignoreNULL = FALSE)

    append_suggestions <- function(chat, num=3) {
      ff <- chat$chat(paste("Suggest",num,"short follow-up questions. Just return the questions as clean enumerated list, do not use bold or italics."))
      qq <- sub("1. ","",strsplit(ff, split='\n[2-9][.] ')[[1]])
      qq <- trimws(qq)
      qq <- paste("  <li class='suggestion submit'>",qq,"</li>")
      qq <- paste("<ul>\n",paste(qq,collapse='\n'),"\n</ul>")
      msg <- list( role = "assistant", content = qq)
      shinychat::chat_append_message("chat", msg, chunk=TRUE, operation="append")
    }
    
    ask_copilot <- function(question, showq=TRUE, suggest=TRUE) {
      if(is.null(chat)) return(NULL)
      if(n_turns() > maxturns) {
        shinyalert::shinyalert(
          title = "Sorry...",
          text = paste("You reached your",maxturns,"token limit. Please upgrade for more tokens."),
          immediate = TRUE,
          showCancelButton = FALSE,
          showConfirmButton = TRUE
        )
        return(NULL)
      }
      if(showq) {
        msg <- list(role="user", content=question)
        shinychat::chat_append_message("chat", msg, chunk=FALSE)
      }
      response <- chat$chat_async(question)
      shinychat::chat_append("chat", response) %...>% {
        if(suggest) {
          append_suggestions(chat, num=2) 
        }
        n_turns(n_turns()+1)
      }
    }

    ## ----------------- tools -------------------------------
    register_tools <- function(chat) {
      if(is.null(chat)) return(NULL)
      chat$register_tool( playbase::ai.tool_get_current_time )
      chat$register_tool( playbase::ai.tool_get_expression )
      #chat$register_tool( playbase::ai.tool_plot_volcano )
    }

    ## ------------ MAIN USER INTERACTION LOOP --------------------
    session$onFlushed( function() {
      dbg("[CopilotModule] onFlushed")
      mesg <- paste("👋 Hi, I'm **BigOmics Copilot**! Ask me about your data")
      shinychat::chat_append("chat", mesg)
    }, once = TRUE)
        
    observeEvent( input$chat_user_input, {
      req(!is.null(chat))
      ask_copilot(input$chat_user_input, showq=FALSE, suggest=input$followup)
    })
        
    ## ---------------------- examples -----------------------------
    observeEvent(input$ask_describe, {
      req(chat)
      ask_copilot("Describe my experiment", suggest=input$followup)
    })
    
    observeEvent(input$ask_findings, {
      req(chat)
      ask_copilot("Summarize main findings of this experiment", suggest=input$followup)
    })
    
    observeEvent(input$ask_pathways, {
      req(chat)
      ask_copilot("List the top enriched pathways and elaborate how they make sense in light of the experiment",
        suggest=input$followup)
    })

    observeEvent(input$show_biomarkers, {
      req(chat)
      ask_copilot("Show the top 3-5 candidate biomarkers (up and down) for this experiment", suggest=input$followup)
    })

    observeEvent(input$find_references, {
      req(chat)
      ask_copilot("Find literature references that are relevant to my experiment", suggest=input$followup)
    })
    
    observeEvent(input$get_expression, {
      req(chat)
      ask_copilot("Get the expression values of MTOR", suggest=input$followup)
    })

    observeEvent(input$plot_volcano, {
      req(chat)
      ask_copilot("Using the tools, create a volcano plot for contrast=1. Return one-line code",
        suggest=FALSE)
    })

    ## ----------------- plot output -------------------------------
    output$plot <- renderPlot({
                  
      if(is.null(chat)) {
        plot.new()
        return()
      }

      if(is.null(chat$last_turn())) return(NULL)
      
      ## parse for R code in last response
      last_response <- chat$last_turn()@text
      m <- regexpr("```[rR].+?```", last_response)
      rcode <- regmatches(last_response, m)

      plotcode <- NULL
      if(length(rcode)) {

        dbg("[CopilotModule:output$plot] WARNING: r code detected!")
        dbg("[CopilotModule:output$plot] rcode = ", rcode)
        
        ## parse for valid plotting code        
        if(grepl("barplot",rcode)) {
          m2 <- regexpr("barplot\\(.+\\)", rcode)
          plotcode <- regmatches(rcode, m2)
        }
        if(grepl("playbase::pgx.Volcano",rcode)) {
          m2 <- regexpr("playbase::pgx.Volcano\\(.+?\\)", rcode)
          plotcode <- regmatches(rcode, m2)
        }
        if(length(plotcode) && nchar(plotcode) > 0) {
          dbg("WARNING: running PLOTCODE = ", plotcode)
          ## super dangerous!!!
          eval(parse(text=plotcode))  
          ##plot(mtcars[,1:3])
        }
      }
      
    })
    
  })
}
