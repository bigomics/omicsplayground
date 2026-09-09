##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
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

  OmicsBoardUI(
    ns = ns,
    title = "Summary",
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
    OmicsBoard(session, pgx, title = "Summary", infotext = "Summary")

    ## pgx$wgcna$report$bullets/$infographic (the old ai-old-wgcnareport.R
    ## schema) is never populated by the current pipeline. The v1 AI-report
    ## schema's "combined" slot -- labelled "Summary" in AI Studio's own
    ## infographic tab (.AI_INFOGRAPHIC_LABELS) -- is the dataset-wide
    ## equivalent: pgx$ai$combined$report (a markdown write-up) and
    ## pgx$ai$combined$infographic (raw PNG bytes + a status/error, from
    ## playbase::pgx.update_infographics(), ai-infographics.R).
    output$bullet_points <- shiny::renderUI({
      txt <- pgx$ai$combined$report
      if(is.null(txt)) txt <- paste(c(paste("- Bullet", 1:5),""),collapse="\n")
      tagList(
        shiny::HTML(markdown::markdownToHTML(txt, fragment.only=TRUE))
      )
    })

    image.RENDER <- function() {
      img.src <- NULL
      infographic <- pgx$ai$combined$infographic
      has.infographic <- !is.null(infographic) &&
        identical(infographic$status, "done") &&
        length(infographic$bytes) > 0
      ## shiny::validate(shiny::need(has.infographic, "Infographic not available."))
      if(has.infographic) {
        img.src <- tempfile(fileext=".png")
        writeBin(infographic$bytes, img.src)
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

