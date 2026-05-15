##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

InfographicSettings <- function(id) {
  ns <- shiny::NS(id)
  info ="AI reports for the current dataset. These reports are AI-generated and can be inaccurate; please always double-check its responses."
  shiny::div(
    style = "padding: 10px 15px;",
    h3("Infographics"),
    shiny::br(),
    div(info),
    shiny::br(),
    br(),
    br(),
    shiny::actionButton(ns("generate"), "Generate images",
      class = "btn btn-primary", style="margin-bottom: 6px;"),
    br(),
    shiny::checkboxInput(ns("force"), "force", FALSE)
  )
}

InfographicUI <- function(id) {
  ns <- shiny::NS(id)
  ui <- bslib::navset_tab(
    id = ns("navset"),
    #bslib::nav_panel(title = "DE", p("DE tab content.")),
    #bslib::nav_panel(title = "Enrichment", p("Enrichment tab content.")),
    bslib::nav_panel(title = "WGCNA",
      bslib::card(
        min_height = 600,
        full_screen = TRUE,
        div(shiny::imageOutput(ns("wgcna"), height="100%", width="100%"),
          height="100%", width="100%", style = "text-align: center;")
      )
    ),
    bslib::nav_panel(title = "moxWGCNA",
      bslib::card(
        min_height = 600,
        full_screen = TRUE,
        div(shiny::imageOutput(ns("wgcna2"), height="100%", width="100%"),
          height="100%", width="100%", style = "text-align: center;")
      )
    ),
    bslib::nav_panel(title = "L1000",
      bslib::card(
        min_height = 600,
        full_screen = TRUE,
        div(shiny::imageOutput(ns("cmap"), height="100%", width="100%"),
          height="100%", width="100%", style = "text-align: center;")
      )
    ),
    bslib::nav_panel(title = "MOFA",
      bslib::card(
        min_height = 600,
        full_screen = TRUE,
        div(shiny::imageOutput(ns("mofa"), height="100%", width="100%"),
          height="100%", width="100%", style = "text-align: center;")
      )
    )
  )    
  return(ui)
}

InfographicServer <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    tmpdir <- tempdir()

    shiny::observeEvent( pgx$X, {
      nav_toggle <- function(x,target) {
        if(is.null(x)) bslib::nav_hide("navset",target)
        if(!is.null(x)) {
          bslib::nav_show("navset",target)
          bslib::nav_select("navset",target)          
        }
      }
      nav_toggle(pgx$gset.meta,"Enrichment")      
      nav_toggle(pgx$mofa,"MOFA")      
      nav_toggle(pgx$gx.meta,"DE")
      nav_toggle(pgx$drugs,"L1000")
      nav_toggle(pgx$wgcna_mox,"moxWGCNA")      
      nav_toggle(pgx$wgcna,"WGCNA")
    })
    
    ## ------------- generate infographics ---------------
    shiny::observeEvent( input$generate, {

      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Please wait. Rendering infographics...", value = 0.3)

      img_model = "google:gemini-3.1-flash-image-preview"
      llm_model = "groq:openai/gpt-oss-120b"
      img_model <- getUserOption(session, "img_model")
      llm_model <- getUserOption(session, "llm_model")      

      dbg("[InfographicServer] llm_model = ", llm_model)
      dbg("[InfographicServer] img_model = ", img_model)      
      
      if (is.null(llm_model) || is.null(img_model)) return(NULL)     
      pgx <- playbase::pgx.update_infographics(
        pgx, force = input$force,
        llm_model=llm_model, img_model=img_model)       
    }) 

    ##--------------- outputs -----------------
    output$wgcna <- renderImage({
      img <- pgx$wgcna$report$infographic
      shiny::validate(need(!is.null(img),"missing WGCNA infographic"))      
      target <- file.path(tmpdir, "infographic-wgcna.png")
      png::writePNG(img, target = target)
      list(src = target, height = "auto", width = "100%")
    }, deleteFile = FALSE)

    output$wgcna2 <- renderImage({
      img <- pgx$wgcna_mox$report$infographic
      shiny::validate(need(!is.null(img),"missing moxWGCNA infographic"))      
      target <- file.path(tmpdir, "infographic-wgcna2.png")
      png::writePNG(img, target = target)
      list(src = target, height = "auto", width = "100%")
    }, deleteFile = FALSE)

    output$cmap <- renderImage({
      img <- pgx$drugs[[1]]$report$infographic
      shiny::validate(need(!is.null(img),"missing cmap infographic"))
      target <- file.path(tmpdir, "infographic-cmap.png")
      png::writePNG(img, target = target)
      list(src = target, height = "auto", width = "100%")
    }, deleteFile = FALSE)

    output$mofa <- renderImage({
      img <- pgx$mofa$report$infographic
      shiny::validate(need(!is.null(img),"missing MOFA infographic"))      
      target <- file.path(tmpdir, "infographic-mofa.png")
      png::writePNG(img, target = target)
      list(src = target, height = "auto", width = "100%")
    }, deleteFile = FALSE)
    
    
  }) ## end of moduleServer
}
