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
      ns("what2show"), "Show:", c("report", "prompt"),
      selected = "report", inline = TRUE
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
    ## shiny::actionButton(
    ##   ns("generate_diagram"), "Regenerate",
    ##   icon = icon("refresh"),
    ##   class = "btn-outline-primary",
    ##   width = "100%"
    ## ),
    shiny::radioButtons(
      ns("diagram_layout"), "Layout:", c("TB", "LR"),
      selected = "TB", inline = TRUE, width = "100%"
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
    shiny::checkboxInput(ns("use_diagram"), "Use diagram", TRUE),
    shiny::selectInput(ns("img_model"), "AI model:", choices = img_models),
    shiny::actionButton(ns("generate_infographic"), "Regenerate", icon = icon("refresh"))
  )

  PlotModuleUI(
    ns("infographic"),
    plotlib = "image",
    title = title,
    label = label,
    #options = options,
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
    if (opt$DEVMODE) {
      shiny::sliderInput(
        ns("topratio"),
        "Top ratio:", 0.6, 0.95, 0.85, 0.05
      )
    },
    if (opt$DEVMODE) {
      shiny::textAreaInput(
        ns("userprompt"),
        label = "User prompt:",
        value = "Be concise. Do not use tables. Use prose. Conclude with a discussion.",
        rows = 5
      )
    },
    ## shiny::actionButton(
    ##   ns("generate_btn"), "Generate!",
    ##   icon = icon("refresh"),
    ##   class = "btn-outline-primary",
    ##   width = "100%"
    ## ),
    shiny::downloadButton(
      ns("downloadPDF"),
      label = "Download",
      class = "btn-outline-primary",
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
    generate_btn <- reactiveVal(0)

    observeEvent( wgcna(), {
      generate_btn(runif(1))
      #shinyjs::disable("downloadPDF") 
    })

    ## observeEvent(input$generate_btn, {
    ##   generate_btn(generate_btn() + 1)
    ##   shinyjs::disable("downloadPDF") 
    ## })

    get_report <- shiny::reactive({
      this_wgcna <- wgcna()
      rpt <- this_wgcna$report
      shinyjs::disable("downloadPDF") 
      validate(need(isTruthy(rpt),"Dataset has no report. Please create with Studio."))
      shinyjs::enable("downloadPDF") 
      return(rpt)
    })
    
    ## ----------------------------------------------------------------------
    ## ------------------------- text module --------------------------------
    ## ----------------------------------------------------------------------

    output$report_bullets <- shiny::renderUI({
      rpt <- get_report()
      txt <- "Generate highlights..."
      if (!is.null(rpt$bullets) && rpt$bullets != "") txt <- rpt$bullets
      txt <- markdown::markdownToHTML(txt, fragment.only = TRUE)
      txt <- gsub("<p>|</p>", "", txt) ## remove p
      shiny::HTML(txt)
    })

    contents_text <- shiny::reactive({
      rpt <- get_report()
      validate(need(isTruthy(rpt),"Dataset has no report. Please create with Studio."))
      if (input$what2show == "prompt") {
        q <- rpt$report_prompt
        txt <- paste("\n\n***Prompt***\n\n", q, "\n")
      }
      if (input$what2show == "report") {
        txt <- rpt$report
      }
      return(txt)
    })

    text.RENDER <- function() {
      txt <- contents_text()
      if(is.null(txt)) shinyjs::disable("downloadPDF") 
      res <- markdown::markdownToHTML(txt, fragment.only = TRUE)
      out <- shiny::div(class = "gene-info", shiny::HTML(res))
      out
    }

    text.RENDER2 <- function() {
      txt <- contents_text()
      res <- markdown::markdownToHTML(txt, fragment.only = TRUE)
      shiny::div(shiny::HTML(res), class = "gene-info", style = "font-size:16px;")
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

    ## ----------------------------------------------------------------------
    ## --------------------------- download ---------------------------------
    ## ----------------------------------------------------------------------

    output$downloadPDF <- downloadHandler(
      filename = function() {
        "wgcna-report.pdf"
      },
      content = function(file) {
        rpt <- get_report()

        dbg("[WGCNA:output$downloadPDF] 1:")
        
        if (is.null(rpt)) {
          playbase::markdownToPDF("PDF report not ready", file = file)
          return()
        }

        dbg("[WGCNA:output$downloadPDF] 2:")
        
        ## replace with updated report, infographic and diagram
        imgpath <- get_infographic()
        rpt$infographic <- imgpath
        rpt$diagram <- get_diagram()
        this_wgcna <- wgcna()
        this_wgcna$report <- rpt

        dbg("[WGCNA:output$downloadPDF] 3:")
        
        ## compile full report
        full_rpt <- playbase::rpt.compile_wgcna_report(
          this_wgcna,
          report = rpt
        )
        playbase::markdownToPDF(full_rpt, file = file)
      }
    )

    ## ----------------------------------------------------------------------
    ## ---------------------------- diagram ---------------------------------
    ## ----------------------------------------------------------------------

    get_diagram <- reactive({
      rpt <- get_report()
      dg <- rpt$diagram
      return(dg)
    })

    diagram.RENDER <- function() {
      dg <- get_diagram()
      shiny::validate(shiny::need(isTruthy(dg), "Diagram not available"))
      lt <- input$diagram_layout
      if (lt == "LR") dg <- sub("rankdir=TB", "rankdir=LR", dg) ## change layout
      if (lt == "TB") dg <- sub("rankdir=LR", "rankdir=TB", dg) ## change layout
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

    ## ----------------------------------------------------------------------
    ## ---------------------------- infographic -----------------------------
    ## ----------------------------------------------------------------------

    empty_infographic <- function(msg="missing infographic") {
      outfile <- tempfile(fileext = ".png")
      png(outfile, w = 1400, h = 800)
      plot.new()
      text(0.5, 0.35, msg, cex = 2.5)
      dev.off()
      return(outfile)
    }
    
    get_infographic <- reactive({
      rpt <- get_report()
      img <- rpt$infographic
      if(!is.null(img)) {
        imgpath <- tempfile(pattern="infographic-", fileext=".png")
        png::writePNG(img, target=imgpath)
      } else {
        ##shiny::validate(shiny::need(isTruthy(img), "Infographic not available"))
        imgpath <- empty_infographic()
      }
      return(imgpath)
    })

    infographic.RENDER <- function() {
      imgpath <- get_infographic()
      ##shiny::validate(shiny::need(!is.null(imgpath), "Infographic path not available."))
      list(
        src = imgpath,
        width = "100%",
        height = "100%",
        alt = "WGCNA Infographic"
      )
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
