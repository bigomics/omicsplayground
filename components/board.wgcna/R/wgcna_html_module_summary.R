##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

wgcna_html_module_summary_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height,
    width) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::radioButtons(
      ns("docstyle"), "Ai summary:", c("short","long"), inline=TRUE
    ),
    shiny::checkboxInput(
      ns("show_prompt"), "Show prompt", FALSE      
    ),
    shiny::actionButton(
      ns("refresh_btn"), "Generate!",
      icon = icon("refresh"),
      class = "btn-outline-primary"
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

wgcna_html_module_summary_server <- function(id,
                                             wgcna,
                                             multi = FALSE,
                                             r_annot = reactive(NULL),
                                             r_module,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    get_summary <- function() {
      wgcna <- wgcna()
      module <- r_module()
      shiny::req(wgcna)      
      res <- "Summary not available"

      if(multi) {
        summaries <- lapply(wgcna,function(w) names(w$summary))
        has.summary <- any(sapply(summaries, function(s) module %in% s))
        if(has.summary) {
          k <- which(sapply(summaries, function(s) module %in% s))
          res <- wgcna[[k]]$summary[[module]]
        }
      } else {
        if("summary" %in% names(wgcna) && module %in% names(wgcna$summary)) {
          res <- wgcna$summary[[module]]
        }
      }
      return(res)
    }

    ##
    ##
    ##
    
    btn_count <- reactiveVal(0)

    observe({
      wgcna <- wgcna()
      module <- r_module()
      btn_count(runif(1))
    })

    observeEvent(input$refresh_btn, {
      btn_count( btn_count() + 1)
    })
    
    contents_text <- shiny::eventReactive({
      btn_count()
    } ,{
      wgcna <- wgcna()
      module <- r_module()        
      ai_model <- getUserOption(session,'llm_model')

      if( btn_count() > 1 && ai_model!='' ) {
        
        docstyle <- switch( input$docstyle,
          "short" = "short summary",
          "long" = "long detailed scientific discussion"          
        )

        annot <- r_annot()
        if(is.null(annot) && !is.null(wgcna$annot)) {
          annot <- wgcna$annot
        }
        
        res <- playbase::wgcna.describeModules(
          wgcna,
          modules = module,
          model = ai_model,
          multi = multi,
          annot = annot,
          docstyle = docstyle,
          numpar = 2,
          experiment = wgcna$experiment,
          verbose=0
        )
        txt <- res$answers[[1]]
        if(input$show_prompt) {
          q <- res$questions[[1]]
          txt <- paste(txt, "\n\n**Prompt**\n\n",q,"\n")
        }
      } else {
        txt <- get_summary()
      }
      txt <- gsub(module,paste0("**",module,"**"),txt)
      txt <- markdown::markdownToHTML(txt, fragment.only=TRUE)
      txt <- paste0("<b>",module," module</b><br><br>", txt)
      return(txt)
    },
    ignoreNULL = FALSE,
    ignoreInit = FALSE
    )
    
    
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


  })
}
