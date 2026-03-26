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
    ##shiny::radioButtons(ns("format"),"format",c("html","PDF"), inline=TRUE)
    shiny::checkboxInput(ns("as_html"),"as html",TRUE)
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

drugconnectivity_report_infographic_ui <- function(
  id,
  label = "",
  title = "",
  info.text = "",
  caption = "",
  height,
  width) {
  ns <- shiny::NS(id)
  
  img_models <- playbase::ai.get_image_models()   
  options <- shiny::tagList(
    shiny::selectInput(ns("img_model"),"AI model:",choices=img_models),
    shiny::actionButton(ns("generate_infographic"),"Generate")
  )

  PlotModuleUI(
    ns("infographic"),
    plotlib = "image",
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


drugconnectivity_report_inputs <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::textAreaInput(
      ns("ai_prompt"),
      label = "User prompt:",
      value = "Be concise. Do not use tables. Use prose. Include a discussion.",
      rows = 5
    ),
    shiny::actionButton(
      ns("ai_generate"), "Generate!",
      icon = icon("refresh"),
      class = "btn-outline-primary",
      width = '100%'
    ),
    shiny::downloadButton(
      ns("downloadPDF"),
      label = "Download PDF",
      class = "btn-outline-primary",      
      width = '100%')
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
                                           rdb,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    btn_count <- reactiveVal(0)
    
    observeEvent( drugs(), {
      btn_count(runif(1))
      infographic_info(" ")      
    })

    observeEvent(input$ai_generate, {
      btn_count( btn_count() + 1)
    })
    
    get_report <- shiny::eventReactive({
      btn_count()
    } ,{
      
      ##llm_model <- "groq:openai/gpt-oss-20b"
      llm_model <- getUserOption(session,'llm_model')
      if(btn_count() < 1 || llm_model == '') {
        return(NULL)
      }
      
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "creating AI report...", value = 0)

      pgx$drugs <- drugs()
      db <- rdb()
      rpt <- playbase::ai.create_report_drug_connectivity(
        pgx, model=llm_model, db=db,
        user.prompt = input$ai_prompt)
        
      return(rpt)
    },
    ignoreNULL = FALSE,
    ignoreInit = FALSE
    )

    output$downloadPDF <- shiny::downloadHandler(
      filename = function() {
        "drugcmap-report.pdf"
      },
      content = function(file) {
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "exporting to PDF...", value = 0.33)
        rpt <- get_report()
        rpt$infographic <- infographic_path()
        db <- rdb()
        full_rpt <- playbase::rpt.compile_drugconnectivity_report(
          obj=NULL, which.db=db, report=rpt,
          pgx=pgx, model=NULL, hlevel=2, shift=TRUE)
        playbase::markdownToPDF(full_rpt, file) 
      }
    )    
    
    ##----------------------------------------------------------------------
    ##------------------------- summary module -----------------------------
    ##----------------------------------------------------------------------

    output$report_bullets <- shiny::renderUI({
      rpt <- get_report()
      ##shiny::req(rpt$bullets)
      txt <- "Summarize this report..."
      if(!is.null(rpt$bullets) && rpt$bullets!="") txt <- rpt$bullets
      tagList(
        shiny::HTML(markdown::markdownToHTML(txt, fragment.only=TRUE))
      )
    })
        
    text.RENDER <- function() {
      rpt <- get_report()
      txt <- rpt$report
      shiny::validate(shiny::need(!is.null(txt), "Please enable AI and generate report."))
      bb  <- markdown::markdownToHTML(rpt$bullets, fragment.only=TRUE)
      if(input$as_html) {
        txt <- markdown::markdownToHTML(txt, fragment.only=TRUE)
        txt <- shiny::HTML(txt)
      }
      shiny::tagList(
        ##bs_alert(shiny::HTML(bb), translate=FALSE, closable=FALSE),
        ##shiny::br(),
        shiny::div(txt)
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

    ##----------------------------------------------------------------------
    ##---------------------------- infographic -----------------------------
    ##----------------------------------------------------------------------

    # Store and trigger the generated image path
    infographic_path <- reactiveVal(NULL)

    infographic_info <- function(msg) {
      outfile <- tempfile(fileext = '.png')
      png(outfile,w=600,h=600)
      plot.new()
      text(0.5,0.5, msg, cex=2.5)
      dev.off()
      infographic_path(outfile)
    }
    
    # Create ExtendedTask for background image generation
    infographic_task <- ExtendedTask$new(function(report, model, outfile) {
      future_promise({
        outfile <- try( playbase::ai.cmap_create_infographic(
          report = report,
          model = model,
          filename = outfile,
          aspectRatio = c("4:3","16:9","3:4")[2],          
          add.fallback = TRUE
        ))
        if(inherits(outfile,"try-error")) return(NULL)
        outfile
      })
    }) 
    
    # Trigger the ExtendedTask when button is clicked
    observeEvent({
      #c(input$generate_infographic, get_report())
      c(input$generate_infographic, get_report())
    }, {
      rpt <- get_report()      
      if(is.null(rpt)) return(NULL)
      has.model <- length(input$img_model)>0 && input$img_model[1]!=""
      shiny::validate(shiny::need(has.model, "Please set your GEMINI_API_KEY"))      
      report <- rpt$report
      model <- input$img_model
      infographic_info("starting...")
      outfile <- tempfile(fileext = '.jpg')
      infographic_task$invoke(report, model, outfile)
    })
    
    # Update reactive value when task status changes. This show a kind
    # of progress message
    observeEvent(infographic_task$status(), {
      status <- infographic_task$status()
      dbg("[generate_infographic] task$status = ",status)
      if(status != "success") {
        msg <- " "
        if(status == "initial") msg <- " "
        if(status == "running") msg <- "generating..."
        infographic_info(msg)
      } else {
        task_result <- infographic_task$result()
        infographic_path(task_result)        
      }
    })

    ## # Update reactive value when task completes
    ## observeEvent(infographic_task$result(), {
    ##   task_result <- infographic_task$result()
    ##   infographic_path(task_result)
    ## })

    infographic.RENDER <- function() {
      #shiny::validate(shiny::need(!is.null(infographic_path()), "Infographic not available."))
      shiny::req(infographic_path())
      # Return a list containing the filename
      list(src = infographic_path(),
        width = "100%",
        height = "100%",
        alt = "Drug Connectivity Infographic")
    }
    
    PlotModuleServer(
      "infographic",
      plotlib = "image",
      func = infographic.RENDER,
      pdf.width = 8, pdf.height = 5,
      res = c(75, 100),
      add.watermark = watermark
    )
    
    
  })
}
