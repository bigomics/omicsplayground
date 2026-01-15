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

      ai_model <- getUserOption(session,'llm_model')
      if(btn_count() < 1 || ai_model=='') return(NULL)
              
      wgcna <- wgcna()
      annot <- r_annot()
      if(is.null(annot) && !is.null(wgcna$annot)) {
        annot <- wgcna$annot
      }
      
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "creating AI report...", value = 0)

      ai_model3 <- c("", ai_model, ai_model)
      rpt <- playbase::wgcna.create_report(
        wgcna, ai_model3, annot=annot, multi=multi,
        format="markdown", ntop=100, verbose=1,
        progress = progress) 
      
      return(rpt)
    },
    ignoreNULL = FALSE,
    ignoreInit = FALSE
    )


    ##----------------------------------------------------------------------
    ##----------------------------------------------------------------------
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
    ##----------------------------------------------------------------------
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
      ai_model <- getUserOption(session,'llm_model')        
      if(diag_count() < 1) {
        dg <- rpt$diagram
      } else if(!is.null(ai_model) && ai_model!='') {
        dg <- playbase::wgcna.create_diagram(rpt$report, ai_model, rankdir="LR") 
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

    
  })
}
