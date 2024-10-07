##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##


ExplainPlotUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::actionLink(
    ns("button"),
    shiny::icon("star"),  ## robot, star, bolt
    class = "text-decoration-none module-btn",
    aria_label = "Explain plot"
  )
}

ExplainPlotModule <- function(id,
                              plot_fun,
                              context = "",
                              chat = NULL,
                              plotlib = "ggplot"
                              ) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    has_warned <- FALSE
    
    is_enabled <- function(){
      has.key <- (Sys.getenv("OPENAI_API_KEY") != "")
      has.key
    }
    
    observeEvent( input$button, {
      if(!is_enabled()) {
        shinyalert::shinyalert(
          title = "",
          text = "Sorry. BigOmics Copilot is not enabled for this account",
          size = "xs",
          showCancelButton = FALSE
        )
        ## beepr::beep(sound=1)
        return(NULL)
      }
      if(!has_warned) {
        shinyalert::shinyalert(
          title = "",
          text = "Warning. The following is AI-generated. It may not always be correct. Your plot and its data will be send to ChatGPT. Do you agree?",
          size = "xs",
          showCancelButton = TRUE,
          callbackR = function(x) {
            if(x==TRUE) explain_plot(chat, plot_fun(), model=NULL) 
          }
        )
        has_warned <<- TRUE
      } else {
        explain_plot(chat, plot_fun(), model=NULL) 
      }
      
    })
    
    #' Convert a plot object to a PNG data URI
    #'
    #' @param p The plot object; currently, plotly and ggplot2 are supported. Note
    #'     that plotly requires Python, {reticulate}, and the PyPI packages {plotly}
    #'     and {kaleido}.
    plot_to_img_content <- function(p) {
      tmp <- tempfile(fileext = ".png")
      on.exit(unlink(tmp))
      switch( plotlib,
             "ggplot" =  plot_to_img_content.ggplot(p, tmp),
             "plotly" =  plot_to_img_content.plotly(p, tmp),
             "base" = plot_to_img_content.base(p, tmp)
             )
      elmer::content_image_file(tmp, resize='low')
    }
    
    plot_to_img_content.base <- function(p, tmp) {
      # Save the plot as an image
      png(file = tmp, width = 512, height = 320)
      p()
      dev.off()
    }

    plot_to_img_content.plotly <- function(p, tmp) {
      # Save the plot as an image
      plotly::save_image(p, tmp, width = 512*0.9, height = 320*0.9, scale=1)
    }
    
    plot_to_img_content.ggplot <- function(p, tmp) {
      ggplot2::ggsave(tmp, p, width = 512, height = 320, units = "px", dpi = 72)
    }
    
    explain_plot <- function(chat, p, model, ..., session = getDefaultReactiveDomain()) {

      img_content <- plot_to_img_content(p)
      img_url <- paste("data:image/png;base64,",img_content@data)
      
      system_prompt_str = paste("You are a chatbot integrated in a data visualization dashboard. You will be asked questions about the following plot, its data and how to interprete this plot. Do not answer questions not related to the plot. Only answer questions about biology, genomics, or statistics. Refuse other questions. Be succint. Use 1, maximum 2 paragraphs. Try to make specific observations if you can, but be conservative in drawing firm conclusions and express uncertainty if you can't be confident. Use following context:\n<CONTEXT>\n",context,"\n</CONTEXT>")
      system_prompt_str = paste("You are a chatbot integrated in a data visualization dashboard. You will be asked questions about the following plot, its data and how to interprete this plot. Do not answer questions not related to the plot. Only answer questions about biology, genomics, or statistics. Refuse other questions. Be succint. Use 1, maximum 2 paragraphs. Try to make specific observations if you can. Sometimes suggest a follow up question (but only questions that you are able to answer). Use following context:\n<CONTEXT>\n",context,"\n</CONTEXT>")
      if(is.null(chat)) {
        chat <- elmer::new_chat_openai(system_prompt = system_prompt_str, model="gpt-4o-mini")
      }

      # chat <- chat$clone()
      chat_id <- paste0(id, "_chatid_", sample.int(1e9, 1))      
      chat_nr <- 1
      
      showModal(modalDialog(
        tags$button(
          type="button",
          class="btn-close d-block ms-auto mb-3",
          `data-bs-dismiss`="modal",
          aria_label="Close"
        ),
        div( class = "text-center mb-4",
            style = "margin-bottom: 20px;",
            tags$img(
              src = img_url,
              style = "max-height: 250px;",
              class = "d-block border mx-auto"
            )
            ),
        shinychat::chat_ui(
          ns(chat_id),
          style = "max-height:calc(75vh - 180px); width: min(800px, 100%);",
          fill = FALSE
        ),
        size = "l",
        easyClose = TRUE,
        title = NULL,
        footer = NULL
      ) |> tagAppendAttributes(style = "--bs-modal-margin: 1.75rem;"))
      
            
      onFlushed(function() {
        ## shinychat::chat_append( ns(chat_id), "👋 Hi, I'm **BigO** your Omics Copilot! I'm here to answer questions about this plot")
        stream <- chat$stream_async("Interpret this plot", img_content)
        Sys.sleep(0.1)
        shinychat::chat_append_stream( ns(chat_id), stream)
        beepr::beep(sound=1) 
      })

      observeEvent( input[[paste0(chat_id, "_user_input")]], {
        msg <- input[[paste0(chat_id, "_user_input")]]
        if( chat_nr > 3) {
          shinychat::chat_append_message( ns(chat_id), list(
            role = "assistant",
            content = "Sorry. I am only allowed to answer only 3 questions."
          ))
          #stream <- chat$stream_async("rephrase: 'sorry, you are asking too many questions'")
          #chat_append_stream( ns(chat_id), stream)
          return(NULL)
        } else {
          stream <- chat$stream_async(msg)
          shinychat::chat_append_stream( ns(chat_id), stream)
        }
        chat_nr <<- chat_nr + 1
      })
    }

  })
}
