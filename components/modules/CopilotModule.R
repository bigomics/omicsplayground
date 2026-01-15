library(dotenv)
library(shiny)
library(shinychat)
#library(elmer)
library(bslib)

CopilotUI <- function(id, layout = c("sidebar","fixed")[1]) {
  ns <- shiny::NS(id)

  #MODELS <- c("gpt-4o-mini","gemma2","gemma2:2b","mistral:7b","llama3.2")
  MODELS <- c("groq:openai/gpt-oss-20b", "groq:llama-3.1-8b-instant",
    "openai:gpt-4.1-nano", "llama3.2:3b","mistral:7b")
  
  sections <- c("description", "dataset_info","compute_settings",
    "differential_expression", "geneset_enrichment",
    "drug_similarity","pcsf_report","wgcna_report")

  sysprompt = "You are computational biologist and will be asked questions about the following experiment. Try to make specific observations if you can. Answer brief and succint."
  
  settings_card <- bslib::navset_underline(
    bslib::nav_panel(      
      "Examples",
      br(),
      "Example questions:",
      br(),
      shiny::actionButton(ns("ask_describe"), "Describe my experiment", width='100%', class="xbtn"),
      shiny::actionButton(ns("ask_findings"), "Summarize main findings", width='100%', class="xbtn"),
      shiny::actionButton(ns("ask_pathways"), "What pathways are involved?", width='100%', class="xbtn"),
      shiny::actionButton(ns("show_biomarkers"), "Show top biomarkers", width='100%', class="xbtn"),
      shiny::actionButton(ns("find_references"), "Find references", width='100%', class="xbtn"),
      shiny::actionButton(ns("get_expression"), "Get expression of MTOR", width='100%', class="xbtn"),
      shiny::actionButton(ns("plot_volcano"), "Show volcano plot", width='100%', class="xbtn")
    ),
    bslib::nav_panel(
      "Settings",
      br(),
      shiny::selectInput(ns("model"), "Model:", choices = MODELS, width='100%'),
      shiny::textAreaInput(ns("prompt"), "Prompt:", value=sysprompt, height=120, width='100%'),
      shiny::checkboxInput(ns("followup"),"Suggest follow-up questions", TRUE),
      shiny::checkboxGroupInput(ns("context"),"Context:", choices = sections,        
        selected = sections, inline = TRUE),
      br(),
      actionButton(ns("load"),"Load model")
    ),
    bslib::nav_panel(
      "PGX",
      br(),
      shiny::fileInput(ns('uploadpgx'), 'Please select pgx', accept = ".pgx"),
      br(),
      ##actionButton(ns("loadpgx"),"Load PGX")
    )
  )

  chat_card <- bslib::card(
      class = "border-0",
      fill = FALSE,
      max_height = "800px",
      shinychat::chat_ui(ns("chat"), width="100%", height = "min(100%,770px)",
        fill = FALSE)
  )

  plot_card <- bslib::card(  
      class = "border-0 pt-5",
      fill = TRUE,
      max_height = "800px",
      plotOutput(ns("plot"))
  )

  if(layout == "fixed") {
    ui <- bslib::layout_columns(
      # col_widths = c(2,6,4),
      col_widths = c(3,9),      
      style = "height: min(90%,700)",
      fill = FALSE,
      settings_card,
      chat_card,
      plot_card
    )
  }
    
  if(layout == "sidebar") {
    ui <- bslib::layout_sidebar(
      style = "padding: 0px;",
      sidebar = bslib::sidebar(
        width = 300,
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
        div( chat_card, class = "p-4 pt-5")
      )
    )
  }

  return(ui)
}


#' Note: pgx needs to be reactiveValues
#'
#' 
CopilotServer <- function(id, pgx, input.click, layout="fixed") {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    chat <- NULL
        
    observeEvent( input$uploadpgx, {
      req(input$uploadpgx)
      pgxfile <- input$uploadpgx$name
      loaded_pgx <- playbase::pgx.load( input$uploadpgx$datapath )
      loaded_pgx <- playbase::pgx.initialize(loaded_pgx)
      for (i in 1:length(loaded_pgx)) {
        pgx[[names(loaded_pgx)[i]]] <- loaded_pgx[[i]]
      }
      dbg("[CopilotModule:uploadpgx] pgx$name =  ", pgx$name)
    })

    #' Observe the copilot button and show the modal AI dialog.
    #'
    observeEvent( input.click(), {
      dbg("[CopilotModule] input.click ")
      
      if(is.null(pgx) || is.null(pgx$X)) {
        shinyalert::shinyalert(
          title = "",
          text = "Oops. First load a dataset before using Omics Copilot",
          size = "xs",
          showCancelButton = FALSE
        )
        return(NULL)
      }

      shiny::showModal(modalDialog2(
        title = NULL,
        CopilotUI(id, layout=layout),
        ##footer = NULL,
        size = "xl",
        easyClose = FALSE,
        fade = FALSE
      ))

    })

    ## On new pgx data, reset/create new chatbot
    observeEvent({
      list( input$load, pgx$X )
    }, {

      dbg("[CopilotModule] react on input.load")      
      req( dim(pgx$X), input$model )
      
      content <- playbase::ai.create_report(pgx, sections=input$context, collate=TRUE)
      prompt <- input$prompt
      if(input$followup) {
        prompt <- paste(prompt, "After each answer, suggest 3 short follow-up questions.")
      }
      prompt <- paste(prompt, "Refuse to answer any question that is not about biology and not related to this experiment.")
      prompt <- paste(prompt, "\nThis is the experiment report: <report>", content, "</report>", collapse=" ")

      message("Creating new chat model ", input$model)
      chat <<- playbase::ai.create_ellmer_chat(input$model, system_prompt=prompt)       
      shinychat::chat_append("chat", paste("[Switching to model", input$model, "...]"))

      message("registering tools...")
      register_tools(chat)

    }, ignoreNULL = FALSE)

    ## ----------------- tools -------------------------------
    register_tools <- function(chat) {
      chat$register_tool( playbase::ai.tool_get_current_time )
      chat$register_tool( playbase::ai.tool_get_expression )
      chat$register_tool( playbase::ai.tool_plot_volcano )
    }

    ## ------------ MAIN USER INTERACTION LOOP --------------------
    session$onFlushed( function() {
      dbg("[CopilotModule] onFlushed")
      shinychat::chat_append("chat", paste("👋 Hi, I'm **BigOmics Copilot**! You can ask me about your data"))
    })
    
    ## count number of interactions
    n_turns <- reactiveVal(0)
    
    observeEvent( input$chat_user_input, {
      dbg("[CopilotModule] input.user_input reacted -> sending to chat ")
      req(!is.null(chat))
      response <- chat$chat_async(input$chat_user_input)
      shinychat::chat_append("chat", response) %...>% {
        n_turns( n_turns() + 1 )
      }
    })
        
    ## ---------------------- examples -----------------------------
    ask_copilot <- function(question) {
      msg <- list(role = "user", content = question)
      shinychat::chat_append_message("chat", msg, chunk=FALSE)
      response <- chat$chat_async(msg$content)
      shinychat::chat_append("chat", response) %...>% {
        n_turns(n_turns()+1)
      }
    }

    observeEvent(input$ask_describe, {
      req(chat)
      ask_copilot("Describe my experiment")
    })
    
    observeEvent(input$ask_findings, {
      req(chat)
      ask_copilot("Summarize the main findings")      
    })
    
    observeEvent(input$ask_pathways, {
      req(chat)
      ask_copilot("List the top enriched pathways and elaborate how they make sense in light of the experiment")
    })

    observeEvent(input$show_biomarkers, {
      req(chat)
      ask_copilot("Show the top candidate biomarkers for this experiment")
    })

    observeEvent(input$find_references, {
      req(chat)
      ask_copilot("Find literatur references that are relevant to my experiment")
    })
    
    observeEvent(input$get_expression, {
      req(chat)
      ask_copilot("Get the expression values of MTOR")
    })

    observeEvent(input$plot_volcano, {
      req(chat)
      ask_copilot("Using the tools, create a volcano plot for contrast=1. Return one-line code")
    })

    ## ----------------- plot output -------------------------------
    output$plot <- renderPlot({
      
      dbg("[CopilotModule:output$plot] ***PLOT*** n_turns = ", n_turns())
      dbg("[CopilotModule:output$plot] ***PLOT*** is.null(chat) = ", is.null(chat))
      dbg("[CopilotModule:output$plot] names.pgx = ", names(pgx))
      dbg("[CopilotModule:output$plot] contrasts = ", colnames(pgx$contrasts))      
            
      if(is.null(chat)) {
        ##return(NULL)
        plot.new()
        return()
      }
      dbg("class.last_turn = ", class(chat$last_turn()))
      if(is.null(chat$last_turn())) return(NULL)
      
      last_response <- chat$last_turn()@text
      cat("last_response = ", last_response, "\n")

      ## parse for R code
      m <- regexpr("```[rR].+?```", last_response)
      rcode <- regmatches(last_response, m)

      message("rcode = ", rcode)

      plotcode <- NULL
      if(length(rcode)) {

        dbg("rcode = ", rcode)
        
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
          dbg("SUPERWARNING: running PLOTCODE = ", plotcode)
          ## super dangerous!!!
          eval(parse(text=plotcode))  
          ##plot(mtcars[,1:3])
        }
      }
      
    })
    
  })
}
