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
    shiny::actionButton(ns("generate_infographic"),"regenerate")
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
    shiny::actionButton(
      ns("generate_btn"), "Generate!",
      icon = icon("refresh"),
      class = "btn-outline-primary"
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
      
      ##llm_model <- "groq:openai/gpt-oss-20b"
      llm_model <- getUserOption(session,'llm_model')
      if(btn_count() < 1 || llm_model == '') {
        return(NULL)
      }
      
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "creating AI report...", value = 0)

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

    ##----------------------------------------------------------------------
    ##---------------------------- infographic -----------------------------
    ##----------------------------------------------------------------------

    # Store and trigger the generated image path
    infographic_path <- reactiveVal(NULL)

    infographic_info <- function(msg) {
      outfile <- tempfile(fileext = '.png')
      png(outfile,w=450,h=600)
      plot.new()
      text(0.5,0.5, msg, cex=2.5)
      dev.off()
      infographic_path(outfile)
    }
    
    cmap_create_infographic <- function(report, model, filename,
                                        add.fallback = TRUE)
    {
      prompt <- paste("**Instructions**: Create an infographic for the following pharmacological mechanism-of-action analysis report.\n\n**summary**:", report)
      out <- playbase::ai.create_image_gemini(
        prompt,
        model = model,
        #model = "gemini-2.5-flash-image",
        #model = "gemini-3-pro-image-preview",
        api_key = Sys.getenv("GEMINI_API_KEY"),
        format = "file",
        filename = filename,
        aspectRatio = "3:4",
        imageSize = "1K"
      )
      dbg("[cmap_create_infographic] 2: filename = ",filename)
      dbg("[cmap_create_infographic] 2: out = ",out)      
      if(is.null(out) || !file.exists(out)) return(NULL)
      dbg("[cmap_create_infographic] x:")
      return(out)
    }

    # Create ExtendedTask for background image generation
    infographic_task <- ExtendedTask$new(function(report, model, outfile) {
      ##infographic_info("starting...")
      future_promise({
        #outfile <- tempfile(fileext = '.jpg')
        outfile <- try(cmap_create_infographic(
          report = report,
          model = model,
          filename = outfile,
          add.fallback = TRUE
        ))
        ##if(inherits(outfile,"try-error")) return(NULL)
        outfile
      })
    }) ## |> bslib::bind_task_button("generate_infographic")
    
    # Trigger the ExtendedTask when button is clicked
    observeEvent({
      c(input$generate_infographic, get_report())
    }, {
      rpt <- get_report()      
      if(is.null(rpt)) return(NULL)
      has.model <- length(input$img_model)>0 && input$img_model[1]!=""
      shiny::validate(shiny::need(has.model, "No Gemini image model available. Please set your GEMINI_API_KEY"))      
      report <- rpt$summary
      model <- input$img_model
      dbg("start infographic task...")
      infographic_info("starting...")
      outfile <- tempfile(fileext = '.jpg')
      dbg("outfile = ", outfile)
      dbg("model = ", model)
      infographic_task$invoke(report, model, outfile)
    })
    
    # Update reactive value when task status changes. This show a kind
    # of progress message
    observeEvent(infographic_task$status(), {
      status <- infographic_task$status()
      dbg("[observe:infographic_task$status] status = ", status)
      if(status != "success") {
        msg <- "..."
        if(status == "initial") msg <- "..."
        if(status == "running") msg <- "generating image ..."
        infographic_info(msg)
      } else {
        task_result <- infographic_task$result()
        dbg("[drugconnectivity_report_server] 1: task_result = ", task_result)        
      }
    })

    # Update reactive value when task completes
    observeEvent(infographic_task$result(), {
      task_result <- infographic_task$result()
      dbg("[drugconnectivity_report_server] 2: task_result = ", task_result)
      infographic_path(task_result)
    })

    infographic.RENDER <- function() {
      shiny::validate(shiny::need(!is.null(infographic_path()), "Infographic not available."))
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
