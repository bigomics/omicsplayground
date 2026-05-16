##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

if(0) {
  source("~/Playground/playbase/dev/include.R", chdir=TRUE)
  devtools::load_all("~/Playground/playbase")
  pgx <- playbase::pgx.load("~/Playground/omicsplayground/data/example-data-reports.pgx")
  pgx <- playbase::pgx.load("~/Playground/omicsplayground/data/mox-daf2.pgx")
  pgx <- playbase::pgx.load("~/Playground/omicsplayground/data/klidel-multiomics2.pgx")
  names(pgx)
  names(pgx$wgcna)
  names(pgx$mofa)  
  names(pgx$wgcna$report)
  names(pgx$wgcna_mox$report)  
  names(pgx$mofa)
}

AiReportSettings <- function(id) {
  ns <- shiny::NS(id)

  info ="AI reports for the current dataset. These reports are AI-generated and can be inaccurate; please always double-check its responses."
  
  shiny::div(
    style = "padding: 10px 15px;",
    h3("AI reports"),
    shiny::br(),
    shiny::div(info),
    shiny::br(),
    br(),
    br(),
    shiny::actionButton(ns("generate"), "Generate reports",
      class="btn btn-primary", style="margin-bottom: 6px;"),
    br(),
    shiny::checkboxInput(ns("force"), "force", FALSE)
  )
}

  
AiReportUI <- function(id) {
  ns <- shiny::NS(id)

  ui <- bslib::navset_tab(
    id = ns("navset"),
    bslib::nav_panel(title = "Summary",
      shiny::div( shiny::htmlOutput(ns("summary")), style="align-items: center;")
    ),
    bslib::nav_panel(title = "WGCNA",
      shiny::div( shiny::htmlOutput(ns("wgcna")), style="align-items: center;")
    ),
    bslib::nav_panel(title = "moxWGCNA",
      shiny::div( shiny::htmlOutput(ns("wgcna2")), style="align-items: center;")
    ),
    bslib::nav_panel(title = "L1000",
      shiny::div( shiny::htmlOutput(ns("cmap")), style="align-items: center;")
    ),
    bslib::nav_panel(title = "MOFA", 
      shiny::div( shiny::htmlOutput(ns("mofa")), style="align-items: center;")
    ),
    bslib::nav_panel(title = "DE",
      shiny::div( shiny::htmlOutput(ns("de")), style="align-items: center;")
    ),
    bslib::nav_panel(title = "Enrichment",
      shiny::div( shiny::htmlOutput(ns("enrichment")), style="align-items: center;")
    )
  )

  return(ui)
}

AiReportServer <- function(id, pgx, rnav) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    observe(dbg("[AiReportServer] rnav =", rnav()))
    
    pdf_tempdir <- tempdir()
    shiny::addResourcePath('pdf', pdf_tempdir)
    on.exit(removeResourcePath('pdf'))

    wgcna_pdf <- reactiveVal(NULL)
    wgcna2_pdf <- reactiveVal(NULL)    
    cmap_pdf <- reactiveVal(NULL)
    mofa_pdf <- reactiveVal(NULL)
    de_pdf <- reactiveVal(NULL)
    enrichment_pdf <- reactiveVal(NULL)
    summary_pdf <- reactiveVal(NULL)

    generate_btn <- reactiveVal(0)

    update_nav <- function(pgx) {
      nav_toggle <- function(x,target) {
        if(is.null(x)) bslib::nav_hide("navset",target)
        if(!is.null(x)) {
          bslib::nav_show("navset",target)
          bslib::nav_select("navset",target)          
        }
      }
      nav_toggle(pgx$gset.meta,"Enrichment")
      nav_toggle(pgx$gx.meta,"DE")
      nav_toggle(pgx$mofa,"MOFA")      
      nav_toggle(pgx$drugs,"L1000")
      nav_toggle(pgx$wgcna_mox,"moxWGCNA")      
      nav_toggle(pgx$wgcna,"WGCNA")
      nav_toggle(pgx$report,"Summary")
    }
    
    shiny::observeEvent( pgx$X, {
      wgcna_pdf(NULL)
      wgcna2_pdf(NULL)      
      cmap_pdf(NULL)
      mofa_pdf(NULL)
      de_pdf(NULL)
      enrichment_pdf(NULL)
      update_nav(pgx)
      generate_btn(1)
    })
    
    ## ----------------------- outputs -------------------------
    output$wgcna <- renderUI({
      shiny::validate(need(!is.null(wgcna_pdf()),"missing WGCNA report"))      
      tag <- paste('<iframe style="height: calc(100vh - 200px); width: calc(100% - 30px)" src="',wgcna_pdf(),'"></iframe>', sep = "")
      return(HTML(tag))
    })

    output$wgcna2 <- renderUI({
      shiny::validate(need(!is.null(wgcna2_pdf()),"missing moxWGCNA report"))      
      tag <- paste('<iframe style="height: calc(100vh - 200px); width: calc(100% - 30px)" src="',wgcna2_pdf(),'"></iframe>', sep = "")
      return(HTML(tag))
    })

    output$cmap <- renderUI({
      shiny::validate(need(!is.null(cmap_pdf()),"missing L1000 report"))            
      tag <- paste('<iframe style="height: calc(100vh - 200px); width: calc(100% - 30px)" src="',cmap_pdf(),'"></iframe>', sep = "")
      return(HTML(tag))
    })

    output$mofa <- renderUI({
      shiny::validate(need(!is.null(mofa_pdf()),"missing MOFA report"))
      tag <- paste('<iframe style="height: calc(100vh - 200px); width: calc(100% - 30px)" src="',mofa_pdf(),'"></iframe>', sep = "")
      return(HTML(tag))
    })

    output$de <- renderUI({
      shiny::validate(need(!is.null(de_pdf()),"missing DE report"))            
      tag <- paste('<iframe style="height: calc(100vh - 200px); width: calc(100% - 30px)" src="',de_pdf(),'"></iframe>', sep = "")
      return(HTML(tag))
    })

    output$enrichment <- renderUI({
      shiny::validate(need(!is.null(enrichment_pdf()),"missing Enrichment report"))
      tag <- paste('<iframe style="height: calc(100vh - 200px); width: calc(100% - 30px)" src="',enrichment_pdf(),'"></iframe>', sep = "")
      return(HTML(tag))
    })

    output$summary <- renderUI({
      shiny::validate(need(!is.null(summary_pdf()),"missing summary report"))
      tag <- paste('<iframe style="height: calc(100vh - 200px); width: calc(100% - 30px)" src="',summary_pdf(),'"></iframe>', sep = "")
      return(HTML(tag))
    })
    
    ## ---------------------- generate reports ----------------------
    shiny::observeEvent(input$generate, generate_btn(generate_btn()+1))

    shiny::observeEvent({
      list(generate_btn(), pgx$X)
    }, {

      dbg("[AiReportServer] input.nav = ", input$nav)
      dbg("[AiReportServer] input.navset = ", input$navset)
      
      if(input$force) {
        wgcna_pdf(NULL)
        wgcna2_pdf(NULL)      
        cmap_pdf(NULL)
        mofa_pdf(NULL)
        de_pdf(NULL)
        enrichment_pdf(NULL)
        summary_pdf(NULL)
      }

      progress <- shiny::Progress$new()
      on.exit(progress$close())
      
      ## compute reports (if missing)
      llm_model = "groq:openai/gpt-oss-120b"
      img_model = "google:gemini-3.1-flash-image-preview"
      llm_model <- getUserOption(session, "llm_model")
      img_model <- getUserOption(session, "img_model")      
      img_model <- NULL
      dbg("[AiReportServer] generate_btn = ", generate_btn())
      
      if (!is.null(llm_model) && llm_model != "" && generate_btn()>1) {
        dbg("[AiReportServer] updating reports...")
        dbg("[AiReportServer] llm_model = ", llm_model)
        #dbg("[AiReportServer] img_model = ", img_model)
        progress$set(message = "Please wait. Updating reports...", value = 0.3)
        pgx <- playbase::pgx.update_reports(
          pgx, force=input$force, llm_model, img_model=NULL,
          select = c("wgcna","mofa","cmap","summary") )        
      }

      dbg("[AiReportServer] Converting to PDF...")
      progress$set(message = "Converting to PDF...", value = 0.7)
      
      ## fill panels
      wgcna_rpt <- pgx$wgcna$report$report
      if(!is.null(wgcna_rpt) && is.null(wgcna_pdf())) {
        #wgcna_rpt <- playbase::rpt.compile_wgcna_report(pgx, report=NULL)
        pdf_target <- file.path(pdf_tempdir,"report-wgcna.pdf")
        playbase::markdownToPDF(wgcna_rpt, file=pdf_target, quiet=TRUE)
        wgcna_pdf("pdf/report-wgcna.pdf")
      }

      wgcna2_rpt <- pgx$wgcna_mox$report$report
      if(!is.null(wgcna2_rpt) && is.null(wgcna2_pdf())) {
        pdf_target <- file.path(pdf_tempdir,"report-wgcna2.pdf")
        playbase::markdownToPDF(wgcna2_rpt, file=pdf_target, quiet=TRUE)
        wgcna2_pdf("pdf/report-wgcna2.pdf")
      }
      
      cmap_rpt <- pgx$drugs[[1]]$report$report
      if(!is.null(cmap_rpt) && is.null(cmap_pdf())) {
        pdf_target <- file.path(pdf_tempdir,"report-cmap.pdf")
        playbase::markdownToPDF(cmap_rpt, file=pdf_target, quiet=TRUE)
        cmap_pdf("pdf/report-cmap.pdf")
      }

      mofa_rpt <- pgx$mofa$report$report
      if(!is.null(mofa_rpt) && is.null(mofa_pdf())) {
        pdf_target <- file.path(pdf_tempdir,"report-mofa.pdf")
        playbase::markdownToPDF(mofa_rpt, file=pdf_target, quiet=TRUE)
        mofa_pdf("pdf/report-mofa.pdf")
      }     

      summary_rpt <- pgx$report$report
      if(!is.null(summary_rpt) && is.null(summary_pdf())) {
        pdf_target <- file.path(pdf_tempdir,"report-summary.pdf")
        playbase::markdownToPDF(summary_rpt, file=pdf_target, quiet=TRUE)
        summary_pdf("pdf/report-summary.pdf")
      }     

      update_nav(pgx)
    }) 
    
  }) ## end of moduleServer
}
