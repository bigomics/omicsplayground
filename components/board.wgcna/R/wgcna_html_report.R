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
  
  PlotModuleUI(
    ns("text"),
    outputFunc = htmlOutput,
    title = title,
    label = label,
    info.text = info.text,
#    options = options,
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
      class = "btn-outline-primary"
    ),
    shiny::radioButtons(
      ns("diagram_layout"), "Layout:", c("TB","LR"),
      selected = "TB", inline=TRUE
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

  img_models <- rev(playbase::ai.get_image_models())   
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

wgcna_report_inputs <- function(id) {
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

      wgcna <- wgcna()
      llm_model <- getUserOption(session,'llm_model')
      if(btn_count() < 1 || llm_model == '') {
        ## rpt <- wgcna$report
        ## if(!is.null(rpt) && !is.null(rpt$diagram) ) {
        ##   return(wgcna$report)
        ## } 
        return(NULL)
      }
      
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "creating AI report...", value = 0)

      annot <- r_annot()
      if(is.null(annot) && !is.null(wgcna$annot)) {
        annot <- wgcna$annot
      }

      rpt <- playbase::wgcna.create_report(
        wgcna, llm_model, annot=annot, 
        graph = wgcna$graph, multi=multi,
        format="markdown", ntop=100, verbose=1,
        progress = progress) 
      
      return(rpt)
    },
    ignoreNULL = FALSE,
    ignoreInit = FALSE
    )

    ##----------------------------------------------------------------------
    ##------------------------- text module --------------------------------
    ##----------------------------------------------------------------------
    
    contents_text <- shiny::reactive({
      rpt <- get_report()
      ## shiny::validate(shiny::need(!is.null(rpt), "Report not available. Please enable AI and generate."))
      if (is.null(rpt)) {
        txt <- "Report not available. Please enable AI and generate."
        return(txt)
      }
      if(input$what2show == "prompt") {
        q <- rpt$report_prompt
        txt <- paste("\n\n***Prompt***\n\n",q,"\n")
      }
      if(input$what2show == "report") {
        txt <- rpt$report
      }
      markdown::markdownToHTML(txt, fragment.only=TRUE)
    })
    
    text.RENDER <- function() {
      res <- contents_text()
      shiny::div(class="gene-info", shiny::HTML(res))
    }

    text.RENDER2 <- function() {
      res <- contents_text()
      shiny::div( shiny::HTML(res), style="font-size:22px;" )
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
      if(diag_count() < 1) {
        dg <- rpt$diagram
      } else if(!is.null(llm_model) && llm_model!='') {
        dg <- playbase::wgcna.create_diagram(rpt$report, llm_model, rankdir="LR") 
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
          filename = outfile
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
