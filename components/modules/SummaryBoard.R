##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

SummaryBoardInputs <- function(id) {
  ns <- shiny::NS(id)  

  bigdash::tabSettings(
    br(),
    withTooltip(
      actionLink(ns("options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top",
      options = list(container = "body")
    )
  )  
}

SummaryBoardUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  div(
    boardHeader(
      title = "Summary",
      info_link = ns("summary_info")
    ),
    bslib::layout_columns(
      col_widths = 12,
      height = "calc(100vh - 124px)",
      row_heights = list("auto",1),
      div(class="alert alert-primary p-2 mt-3", 
        htmlOutput(ns("bullet_points"))),
      bslib::layout_columns(
        col_widths = 12,
        height = "100%",
        PlotModuleUI(
          ns("image"),
          plotlib = "image",
          title = "Summary",
          label = "",
          options = NULL,
          info.text = "info text",
          caption = "caption",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%"),
          download.fmt = c("png", "pdf", "svg")
        )
      )
    ) ## end layout cols
  )
}

SummaryBoard <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    bigdash.unloadSidebar()

    fullH <- 800
    tabH <- "70vh"
    
    ## ========================================================================
    ## ============================ OBSERVERS =================================
    ## ========================================================================
    summary_info <- "Summary"

    observeEvent(input$summary_info, {
      showModal(
        modalDialog(
          title = tags$strong("Summary Board"),
          summary_info,
          easyClose = TRUE,
          size = "l"
        )
      )
    })

    output$bullet_points <- shiny::renderUI({
      txt <- pgx$wgcna$report$bullets
      if(is.null(txt)) txt <- paste(c(paste("- Bullet", 1:5),""),collapse="\n")
      tagList(
        shiny::HTML(markdown::markdownToHTML(txt, fragment.only=TRUE))
      )
    })

    image.RENDER <- function() {
      img.src <- NULL
      ## shiny::validate(shiny::need(!is.null(img.src), "Infographic not available."))
      img <- pgx$wgcna$report$infographic
      dim(img)
      if(!is.null(img)) {
        img.src <- tempfile(fileext=".png")
        png::writePNG(img, target=img.src)
      }
      list(
        src = img.src,
        width = "100%",
        height = "100%",
        alt = "Infographic"
      )
    }
    
    PlotModuleServer(
      "image",
      plotlib = "image",
      func = image.RENDER,
      pdf.width = 10, pdf.height = 5,
      res = c(75, 100),
      add.watermark = FALSE
    )
    

  })
}

