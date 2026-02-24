##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Mechanism-of-action plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
drugconnectivity_report_summary_ui <- function(
  id,
  label = "",
  title = "",
  info.text = "",
  caption = "",
  height,
  width
) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::radioButtons(ns("format"),"format",c("html","PDF"), inline=TRUE)
  )
  
  PlotModuleUI(
    ns("summary"),
    outputFunc = htmlOutput,
    title = title,
    label = label,
    info.text = info.text,
    options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

drugconnectivity_report_bullets_ui <- function(id) {
  ns <- shiny::NS(id)
  htmlOutput(ns("report_bullets"))
}

drugconnectivity_report_diagram_ui <- function(
  id,
  label = "",
  title = "",
  info.text = "",
  caption = "",
  height,
  width) {
  ns <- shiny::NS(id)
  
  PlotModuleUI(
    ns("diagram"),
    #outputFunc = imageOutput,
    title = title,
    label = label,
    info.text = info.text,
    ##options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}


drugconnectivity_report_inputs <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::actionButton(
      ns("generate_btn"), "Generate!",
      icon = icon("refresh"),
      class = "btn-outline-primary"
    ),
    shiny::radioButtons(
      ns("what2show"), "Show:", c("report","prompt"),
      selected = "report", inline=TRUE
    )    
  )
}

#' Mechanism of action plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
drugconnectivity_report_server <- function(id,
                                           pgx,
                                           drugs,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    btn_count <- reactiveVal(0)
    
    observeEvent( drugs(), {
      btn_count(runif(1))
    })

    observeEvent(input$generate_btn, {
      btn_count( btn_count() + 1)
    })
    
    get_report <- shiny::eventReactive({
      btn_count()
    } ,{
      
      llm_model <- getUserOption(session,'llm_model')
      if(btn_count() < 1 || llm_model == '') {
        return(NULL)
      }
      
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "creating AI report...", value = 0)

      llm_model <- "groq:openai/gpt-oss-20b"
      pgx$drugs <- drugs()
      rpt <- playbase::ai.summarize_drug_connectivity(pgx, ct=1, db=1, model=llm_model)
      
      return(rpt)
    },
    ignoreNULL = FALSE,
    ignoreInit = FALSE
    )

    ##----------------------------------------------------------------------
    ##------------------------- summary module --------------------------------
    ##----------------------------------------------------------------------

    output$report_bullets <- shiny::renderUI({
      rpt <- get_report()
      ##shiny::req(rpt$bullets)
      txt <- "Summarize this report"
      if(!is.null(rpt$bullets) && rpt$bullets!="") txt <- rpt$bullets
      tagList(
        shiny::HTML(markdown::markdownToHTML(txt, fragment.only=TRUE))
      )
    })
        
    text.RENDER <- function() {
      rpt <- get_report()
      txt <- rpt$summary
      shiny::validate(shiny::need(!is.null(txt), "Please enable AI and generate report."))
      res <- markdown::markdownToHTML(rpt$summary, fragment.only=TRUE)
      bb  <- markdown::markdownToHTML(rpt$bullets, fragment.only=TRUE)
      shiny::tagList(
        ##bs_alert(shiny::HTML(bb), translate=FALSE, closable=FALSE),
        ##shiny::br(),
        shiny::div(shiny::HTML(res))
      )
    }

    PlotModuleServer(
      "summary",
      plotlib = "generic",
      plotlib2 = "generic",
      func = text.RENDER,
      func2 = text.RENDER,
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI,
      pdf.width = 8, pdf.height = 5,
      res = c(75, 100),
      add.watermark = watermark
    )
    
    
  })
}
