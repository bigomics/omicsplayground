##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

wgcna_html_report_ui <- function(
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
    shiny::radioButtons(
      ns("what2show"), "Show:", c("report","prompt"),
      selected = "report", inline=TRUE
    )
  )
  
  PlotModuleUI(
    ns("text"),
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

wgcna_report_diagram_ui <- function(
  id,
  title = "",
  info.text = "",
  caption = "",
  label = "",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::actionButton(
      ns("generate_diagram"), "Regenerate",
      icon = icon("refresh"),
      class = "btn-outline-primary",
      width = "100%"
    ),
    shiny::radioButtons(
      ns("diagram_layout"), "Layout:", c("TB","LR"),
      selected = "TB", inline=TRUE, width = "100%"
    )        
  )

  PlotModuleUI(
    ns("diagram"),
    plotlib = "svgPanZoom",    
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

wgcna_report_infographic_ui <- function(
  id,
  title = "",
  info.text = "",
  caption = "",
  label = "",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  img_models <- playbase::ai.get_image_models(opt$IMAGE_MODELS)   
  options <- shiny::tagList(
    shiny::selectInput(ns("img_model"),"AI model:",choices=img_models),
    shiny::actionButton(ns("generate_infographic"),"regenerate")
  )

  PlotModuleUI(
    ns("infographic"),
    #outputFunc = uiOutput,
    #plotlib = "generic",
    plotlib = "image",
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width
  )
}

wgcna_report_bullets_ui <- function(id) {
  ns <- shiny::NS(id)
  htmlOutput(ns("report_bullets"))
}

wgcna_report_inputs <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    if(opt$DEVMODE) {
      shiny::sliderInput(
        ns("topratio"),
        "Top ratio:", 0.6, 0.95, 0.85, 0.05
      )
    },
    if(opt$DEVMODE) {
      shiny::textAreaInput(
        ns("userprompt"),
        label = "User prompt:",
        value = "Be concise. Do not use tables. Use prose. Conclude with a discussion.",
        rows = 5
      )
    },
    shiny::actionButton(
      ns("generate_btn"), "Generate!",
      icon = icon("refresh"),
      class = "btn-outline-primary",
      width = "100%"
    ),
    shiny::downloadButton(
      ns("downloadPDF"),
      label = "Download",
      style = "width: 100%;")
  )
}


wgcna_html_report_server <- function(id,
                                     wgcna,
                                     multi = TRUE,
                                     r_annot = reactive(NULL),
                                     watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    btn_count <- reactiveVal(0)

    observe({
      wgcna <- wgcna()
      btn_count(runif(1))
    })

    observeEvent(input$generate_btn, {
      btn_count( btn_count() + 1)
    })
    
    get_report <- shiny::eventReactive({
      btn_count()
    } ,{

      llm_model <- getUserOption(session,'llm_model')
      if(is.null(llm_model) || llm_model == '') return(NULL)
      
      this_wgcna <- wgcna()
      if(btn_count() < 1) {
        rpt <- this_wgcna$report
        if(is.null(rpt)) return(NULL)
        dbg("*** found report in wgcna object")
        dbg("*** names.rpt = ", names(rpt))
        ##return(rpt)
        return(NULL)
      }
      
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "creating AI report...", value = 0)

      annot <- r_annot()
      if(is.null(annot) && !is.null(this_wgcna$annot)) {
        annot <- this_wgcna$annot
      }

      userprompt <- ifelse(is.null(input$userprompt),"",input$userprompt)
      topratio <- ifelse(is.null(input$topratio),0.85,input$topratio)
      
      rpt <- playbase::wgcna.create_report(
        this_wgcna,
        ai_model = llm_model,
        annot = annot, 
        graph = this_wgcna$graph,
        multi = multi,
        userprompt = userprompt,
        ntop = 100,
        topratio = topratio,
        psig = 0.05,
        format = "markdown",
        verbose = 1,
        progress = progress) 
      
      return(rpt)
    },
    ignoreNULL = FALSE,
    ignoreInit = FALSE
    )

    output$downloadPDF <- downloadHandler(
      filename = function() {"wgcna-report.pdf"},
      content = function(file) {
        rpt <- get_report()
        if (is.null(rpt)) {
          playbase::markdownToPDF("PDF report not ready", file=file)
          return()
        }

        img <- infographic_path()
        if(img=='') img <- NULL
        if(!is.null(img) && file.exists(img)) {
          dbg("using created infographic for downloaded PDF report")
        } else {
          dbg("Warning: missing infographic for downloaded PDF report")          
        }
        rpt$infographic <- img
        this_wgcna <- wgcna()
        this_wgcna$report <- rpt
        
        full_rpt <- playbase::rpt.compile_wgcna_report(
          this_wgcna, report = rpt)
        playbase::markdownToPDF(full_rpt, file=file) 
      }
    )
    
    ##----------------------------------------------------------------------
    ##------------------------- text module --------------------------------
    ##----------------------------------------------------------------------

    output$report_bullets <- shiny::renderUI({
      rpt <- get_report()
      txt <- "Generate report highlights..."
      if(!is.null(rpt$bullets) && rpt$bullets!="") txt <- rpt$bullets
      txt <- markdown::markdownToHTML(txt, fragment.only=TRUE)
      txt <- gsub("<p>|</p>","",txt) ## remove p
      shiny::HTML(txt)
    })
    
    contents_text <- shiny::reactive({
      rpt <- get_report()
      if (is.null(rpt)) return(NULL)
      if(input$what2show == "prompt") {
        q <- rpt$report_prompt
        txt <- paste("\n\n***Prompt***\n\n",q,"\n")
      }
      if(input$what2show == "report") {
        txt <- rpt$report
      }
      return(txt)
    })
    
    text.RENDER <- function() {
      txt <- contents_text()
      shiny::validate(shiny::need(!is.null(txt), "Please enable AI and generate report."))
      res <- markdown::markdownToHTML(txt, fragment.only=TRUE)
      out <- shiny::div(class="gene-info", shiny::HTML(res))
      out
    }

    text.RENDER2 <- function() {
      txt <- contents_text()
      shiny::validate(shiny::need(!is.null(txt), "Please enable AI and generate report."))
      res <- markdown::markdownToHTML(txt, fragment.only=TRUE)
      shiny::div(shiny::HTML(res), class="gene-info", style="font-size:16px;")
    }

    PlotModuleServer(
      "text",
      plotlib = "generic",
      plotlib2 = "generic",
      func = text.RENDER,
      func2 = text.RENDER2,
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI,
      pdf.width = 8, pdf.height = 5,
      res = c(75, 100),
      add.watermark = watermark
    )
    
    ##----------------------------------------------------------------------
    ##---------------------------- diagram ---------------------------------
    ##----------------------------------------------------------------------

    diag_count <- reactiveVal(0)
    
    observe({
      wgcna <- wgcna()
      diag_count(runif(1))
    })

    observeEvent(input$generate_diagram, {
      diag_count( diag_count() + 1)
    })
    
    get_diagram <- reactive({
      rpt <- get_report()
      shiny::validate(shiny::need(!is.null(rpt), "Diagram not available"))
      llm_model <- getUserOption(session,'llm_model')        
      this_wgcna <- wgcna()
      graph <- this_wgcna$graph
      if(diag_count() < 1) {
        dg <- rpt$diagram
      } else if(!is.null(llm_model) && llm_model!='') {
        dg <- playbase::wgcna.create_diagram(
          rpt$report, llm_model, graph = graph,
          rankdir="LR") 
      } else {
        dg <- NULL
      }
      dg
    })

    diagram.RENDER <- function() {
      dg <- get_diagram()
      lt <- input$diagram_layout
      if(lt=="LR") dg <- sub("rankdir=TB","rankdir=LR",dg)  ## change layout
      if(lt=="TB") dg <- sub("rankdir=LR","rankdir=TB",dg)  ## change layout      
      DiagrammeR::grViz(dg)
    }

    diagram.SVG <- function() {
      dg <- diagram.RENDER() 
      img.svg <- DiagrammeRsvg::export_svg(dg)
      pz <- svgPanZoom::svgPanZoom(
        img.svg,
        controlIconsEnabled = TRUE,
        zoomScaleSensitivity = 0.4,
        minZoom = 1,
        maxZoom = 5,
        viewBox = FALSE
        )
      return(pz)
    }
    
    PlotModuleServer(
      "diagram",
      plotlib = "svgPanZoom",      
      func = diagram.SVG,
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
      png(outfile,w=1600,h=800)
      plot.new()
      text(0.5,0.5, msg, cex=2.5)
      dev.off()
      infographic_path(outfile)
    }
    
    # Create ExtendedTask for background image generation
    infographic_task <- ExtendedTask$new(function(report, diagram, model) {
      infographic_info("starting...")
      future_promise({
        outfile <- tempfile(fileext = '.jpg')
        outfile <- try(playbase::wgcna.create_infographic(
          report = report,
          diagram = diagram,
          model = model,
          filename = outfile,
          add.fallback = TRUE
        ))
        if(inherits(outfile,"try-error")) return(NULL)
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
      report <- rpt$report
      diagram <- rpt$diagram
      model <- input$img_model
      dbg("start infographic task...")
      infographic_task$invoke(report, diagram, model)
    })
    
    # Update reactive value when task status changes. This show a kind
    # of progress message
    observeEvent(infographic_task$status(), {
      status <- infographic_task$status()
      dbg("[observe:infographic_task$status] status = ", status)
      if(status != "success") {
        msg <- "..."
        if(status == "initial") msg <- "waiting for diagram..."
        if(status == "running") msg <- "generating image ..."
        infographic_info(msg)
      }
    })

    # Update reactive value when task completes
    observeEvent(infographic_task$result(), {
      result <- infographic_task$result()
      infographic_path(result)
    })

    infographic.RENDER <- function() {
      shiny::validate(shiny::need(!is.null(infographic_path()), "Infographic not available."))
      # Return a list containing the filename
      list(src = infographic_path(),
        width = "100%",
        height = "100%",
        alt = "WGCNA Infographic")
    }
    
    PlotModuleServer(
      "infographic",
      plotlib = "image",
      #plotlib = "generic",
      func = infographic.RENDER,
      pdf.width = 8, pdf.height = 5,
      res = c(75, 100),
      add.watermark = watermark
    )
    
    
  })
}
