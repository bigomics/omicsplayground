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
    shiny::actionButton(ns("clear"), "Clear",
      class="btn btn-warning-outline", style="margin-bottom: 6px;"),
    br(),
    shiny::checkboxInput(ns("force"), "Force regenerate", FALSE),
    br(),
    shiny::fileInput(ns("logofile"), "Replace logo:", accept=c(".png",".jpg"))
  )
}

  
AiReportUI <- function(id) {
  ns <- shiny::NS(id)

  preview_ui <- bslib::navset_tab(
    id = ns("navset_preview"),
    bslib::nav_panel(title = "Summary",
      shiny::div( shiny::htmlOutput(ns("summary")) %>% bigLoaders::useSpinner(),
        style="align-items: center;")
    ),
    bslib::nav_panel(title = "WGCNA",
      shiny::div( shiny::htmlOutput(ns("wgcna"))%>% bigLoaders::useSpinner(),
        style="align-items: center;")      
    ),
    bslib::nav_panel(title = "moxWGCNA",
      shiny::div( shiny::htmlOutput(ns("wgcna2")) %>% bigLoaders::useSpinner(),
        style="align-items: center;")
    ),
    bslib::nav_panel(title = "L1000",
      shiny::div( shiny::htmlOutput(ns("cmap")) %>% bigLoaders::useSpinner(),
        style="align-items: center;")
    ),
    bslib::nav_panel(title = "MOFA",
      shiny::div( shiny::htmlOutput(ns("mofa")) %>% bigLoaders::useSpinner(),
        style="align-items: center;")
    ),
    bslib::nav_panel(title = "DE",
      shiny::div( shiny::htmlOutput(ns("de")) %>% bigLoaders::useSpinner(),
        style="align-items: center;")
    ),    
    bslib::nav_panel(title = "Enrichment",
      shiny::div( shiny::htmlOutput(ns("enrichment")) %>% bigLoaders::useSpinner(),
        style="align-items: center;")
    )
  )

  text_area <- function(id) {
    div( shiny::textAreaInput(ns(id), NULL, value = "",
      height = "calc(100vh - 160px)", width = "100%"),
      id="reportTA", style='font-size: 1em;')
  }
  
  edit_ui <- bslib::navset_tab(
    id = ns("navset_edit"),
    bslib::nav_panel(title = "Summary", text_area("edit_summary")),
    bslib::nav_panel(title = "WGCNA", text_area("edit_wgcna")),
    bslib::nav_panel(title = "moxWGCNA", text_area("edit_wgcna2")),
    bslib::nav_panel(title = "L1000", text_area("edit_cmap")),
    bslib::nav_panel(title = "MOFA", text_area("edit_mofa")),
    bslib::nav_panel(title = "DE", text_area("edit_de")),
    bslib::nav_panel(title = "Enrichment", text_area("edit_enrichment"))
  )

  ui <- bslib::navset_hidden(
    id = ns("navset"),
    bslib::nav_panel(title = "Preview", preview_ui),    
    bslib::nav_panel(title = "Edit", edit_ui),
    header = tagList(
      bslib::nav_spacer(),
      div( bslib::input_switch(ns("edit_mode"), "edit", FALSE),
        style = "position: absolute; width: 100px; right: 0px; font-size: 1.1em;")
      )
  )
  return(ui)
}


AiReportServer <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    
    pdf_tempdir <- tempdir()
    opg.logo = file.path(OPG,"components/app/R/www/bigomics-logo-small.png")
    logopath = file.path(pdf_tempdir,"logo.png")
    file.copy(opg.logo, logopath)
    
    de_pdf <- file.path(pdf_tempdir,"report-de.pdf")
    enrichment_pdf <- file.path(pdf_tempdir,"report-enrichment.pdf")      
    mofa_pdf <- file.path(pdf_tempdir,"report-mofa.pdf")
    wgcna_pdf <- file.path(pdf_tempdir,"report-wgcna.pdf")
    wgcna2_pdf <- file.path(pdf_tempdir,"report-wgcna2.pdf")      
    cmap_pdf <- file.path(pdf_tempdir,"report-cmap.pdf")
    summary_pdf <- file.path(pdf_tempdir,"report-summary.pdf")      

    shiny::addResourcePath('pdf', pdf_tempdir)
    on.exit({
      unlink(de_pdf)
      unlink(enrichment_pdf)
      unlink(mofa_pdf)
      unlink(wgcna_pdf)
      unlink(wgcna2_pdf)
      unlink(cmap_pdf)
      unlink(summary_pdf)      
      removeResourcePath('pdf')
    })
    
    update_nav <- function(pgx) {
      nav_toggle <- function(x,target) {
        if(is.null(x)) {
          bslib::nav_hide("navset_preview",target)
          bslib::nav_hide("navset_edit",target)          
        }
        if(!is.null(x)) {
          bslib::nav_show("navset_preview",target)
          bslib::nav_select("navset_preview",target)          
          bslib::nav_show("navset_edit",target)
          bslib::nav_select("navset_edit",target)          
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

    has.warned=0
    
    shiny::observeEvent( input$edit_mode, {
      if(input$edit_mode) {
        bslib::nav_select("navset_edit", input$navset_preview)                  
        bslib::nav_select("navset", "Edit")
        ## shinyalert::shinyalert(
        ##   text = HTML("Warning. Edits are currently <b>not saved</b> after logout!"),
        ##   html = TRUE,
        ##   closeOnClickOutside = TRUE,
        ##   timer = ifelse(has.warned, 1000, 0)
        ## )
        has.warned <<- has.warned + 1
      } else {
        bslib::nav_select("navset_preview", input$navset_edit)          
        bslib::nav_select("navset", "Preview")
      }
    })
    
    ## ----------------------- PDF outputs -------------------------
    countLetters <- function(text) {
      az <- c(letters, LETTERS)
      cc <- sapply(az, function(s) sum(gregexpr(s,text)[[1]]>0))
      cc[az]
    }

    char_counts <- c( "wgcna" = 0, "wgcna2" = 0, "mofa" = 0, "cmap" = 0,
      "summary" = 0, "de" = 0, "enrichment" = 0)    
    letter_counts = rep(list(countLetters("")), length(names(char_counts)))
    names(letter_counts) <- names(char_counts)
    
    has.changed <- function(what, now) {
      if(is.null(now)) now <- ""
      ct <- nchar(now)
      changed1 <- (!is.null(now) && ct != char_counts[what])
      char_counts[what] <<- ct
      lt <- countLetters(now)
      changed2 <- !all((lt - letter_counts[[what]]) == 0)
      letter_counts[[what]] <<- lt
      return(changed1 || changed2)
    }

    updatePDF <- function(rpt, rptname, pdfname) {
      if(has.changed(rptname, rpt)) {
        shiny::withProgress(message = "Rendering PDF...", value = 0.7, {
          file <- file.path(pdf_tempdir, pdfname)
          playbase::markdownToPDF(rpt, file=file, logo=logopath, quiet=TRUE)
        })
      }
    }

    pdf.iframe <- function(pdfname) {
      HTML(paste0('<embed style="height: calc(100vh - 200px); width: calc(100% - 20px)" src="pdf/',pdfname,'" type="application/pdf" />'))      
    }
    
    output$summary <- renderUI({
      updatePDF(input$edit_summary, "summary", "report-summary.pdf")
      shiny::validate(need(file.exists(summary_pdf), "missing summary report"))
      return(pdf.iframe("report-summary.pdf"))
    })

    output$wgcna <- renderUI({
      updatePDF(input$edit_wgcna, "wgcna", "report-wgcna.pdf")      
      shiny::validate(need(file.exists(wgcna_pdf),"missing WGCNA report"))      
      return(pdf.iframe("report-wgcna.pdf"))
    })

    output$wgcna2 <- renderUI({
      updatePDF(input$edit_wgcna2, "wgcna2", "report-wgcna2.pdf")            
      shiny::validate(need(file.exists(wgcna2_pdf),"missing moxWGCNA report"))      
      return(pdf.iframe("report-wgcna2.pdf"))
    })

    output$cmap <- renderUI({
      updatePDF(input$edit_cmap, "cmap", "report-cmap.pdf")                  
      shiny::validate(need(file.exists(cmap_pdf),"missing L1000 report"))            
      return(pdf.iframe("report-cmap.pdf"))
    })
    
    output$mofa <- renderUI({
      updatePDF(input$edit_mofa, "mofa", "report-mofa.pdf")                        
      shiny::validate(need(file.exists(mofa_pdf),"missing mofa report"))            
      return(pdf.iframe("report-mofa.pdf"))
    })

    output$de <- renderUI({
      ## updatePDF(input$edit_de, "de", "report-de.pdf")                        
      shiny::validate(need(file.exists(de_pdf),"missing DE report"))            
      return(pdf.iframe("report-de.pdf"))
    })

    output$enrichment <- renderUI({
      ## updatePDF(input$edit_enrichment, "de", "report-enrichment.pdf")      
      shiny::validate(need(file.exists(enrichment_pdf),"missing Enrichment report"))
      return(pdf.iframe("report-enrichment.pdf"))
    })

    
    ##---------------------------------------------------------------
    ## --------- main observer: generate reports --------------------
    ##---------------------------------------------------------------
    
    shiny::observeEvent({
      list(input$generate, pgx$X, pgx$name)
    }, {
      shiny::req(pgx$X, pgx$name)
      ##shiny::req(!is.null(input$generate))
      dbg("[[AiReportServer]] is.null(input$generate)", is.null(input$generate))
      
      ## compute reports (if missing)
      llm_model = "groq:openai/gpt-oss-120b"
      img_model = "google:gemini-3.1-flash-image-preview"
      llm_model <- getUserOption(session, "llm_model")
      img_model <- getUserOption(session, "img_model")      
      img_model <- NULL
      if (llm_model=="") llm_model <- NULL

      has.reports <- playbase::pgx.has_reports(pgx)
      dbg("[AiReportServer] has.reports = ", has.reports)
      
      if (!is.null(llm_model) && llm_model != "") {
        dbg("[AiReportServer] updating reports using llm = ", llm_model)
        
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Please wait. Updating reports...", value = 0.3)

        ##pgx.showSmallModal("Please wait. Updating reports...")
        pgx <- playbase::pgx.update_reports(
          pgx, force = input$force, llm_model, img_model = NULL,
          select = c("wgcna","mofa","cmap","summary") )
        ##shiny::removeModal()

      }

      update_nav(pgx)      
      clear_files()

      rpt_wgcna <- pgx$wgcna$report$report
      rpt_wgcna2 <- pgx$wgcna_mox$report$report      
      rpt_mofa <- pgx$mofa$report$report
      rpt_cmap <- pgx$drugs[[1]]$report$report      
      rpt_summary <- pgx$report$report
      
      updateTextAreaInput(session, "edit_wgcna", value = rpt_wgcna)
      updateTextAreaInput(session, "edit_wgcna2", value = rpt_wgcna2)
      updateTextAreaInput(session, "edit_mofa", value = rpt_mofa)
      updateTextAreaInput(session, "edit_cmap", value = rpt_cmap)
      updateTextAreaInput(session, "edit_summary", value = rpt_summary)

      updateCheckboxInput(session, "force", value = FALSE)
      updateActionButton(session, "generate", label = "Reset reports")      
    })
    
    ##-------------------------------------------------
    ##------------ other observers --------------------
    ##-------------------------------------------------

    observe({
      file <- input$logofile
      req(file)
      ext <- tools::file_ext(file$datapath)
      validate(need(ext %in% c('jpg','png'), "Please upload a png or jpg file"))
      file.copy(file$datapath, logopath, overwrite=TRUE)
    })

    shiny::observeEvent({
      list(pgx$X, input$clear)
    }, {
      clear_files()
      update_nav(pgx)      
      updateCheckboxInput(session, "force", value = FALSE)      
      updateActionButton(session, "generate", label = "Generate reports")
    })

    clear_files <- function() {
      dbg("[clear_files] clearing files!")
      updateTextAreaInput(session, "edit_wgcna", value = "")
      updateTextAreaInput(session, "edit_wgcna2", value = "")
      updateTextAreaInput(session, "edit_mofa", value = "")
      updateTextAreaInput(session, "edit_cmap", value = "")
      updateTextAreaInput(session, "edit_summary", value = "")
      unlink(file.path(pdf_tempdir,"report-wgcna.pdf"))
      unlink(file.path(pdf_tempdir,"report-wgcna2.pdf"))
      unlink(file.path(pdf_tempdir,"report-cmap.pdf"))
      unlink(file.path(pdf_tempdir,"report-mofa.pdf"))
      unlink(file.path(pdf_tempdir,"report-summary.pdf"))            
      char_counts <<- char_counts * 0
      letter_counts <<- lapply(letter_counts, function(x) x*0)
    }

    
  }) ## end of moduleServer
}
