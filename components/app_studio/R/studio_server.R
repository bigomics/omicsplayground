##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' The application server-side logic
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @export
#' Note: pgx needs to be reactiveValues
#'
#'
StudioServer <- function(id, pgx, save_pgx = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    OmicsBoard("board", pgx, title="AI Studio", infotext = NULL) 
    
    ## ------------------- subpanels UI ------------------------
    output$dataset_info <- renderUI({
      HTML("Hello world!")
    })
    
    ## ---------------------- observers ------------------------

    observeEvent(input$show_poster, {
      bslib::nav_select("studiopanel","Poster")
      bslib::nav_select("settings","Poster")      
    })

    observeEvent(input$show_slidedeck, {
      bslib::nav_select("studiopanel","Slide deck")
      bslib::nav_select("settings","Slide deck")      
    })
    
    observeEvent(input$show_reports, {
      bslib::nav_select("studiopanel","Reports")
      bslib::nav_select("settings","Reports")      
    })

    observeEvent(input$show_infographic, {
      bslib::nav_select("studiopanel","Infographic")
      bslib::nav_select("settings","Infographic")      
    })
    
    observeEvent(input$show_podcast, {
      bslib::nav_select("studiopanel","Podcast")
      bslib::nav_select("settings","Podcast")
    })

    observeEvent(input$show_quiz, {
      bslib::nav_select("studiopanel","Quiz")
      bslib::nav_select("settings","Quiz")
    })

    observeEvent(input$show_custom, {
      bslib::nav_select("studiopanel","Custom")
      bslib::nav_select("settings","Custom")
    })
    
    ## ----------------- output -------------------------------
    VisReportServer("poster", pgx, output_format="poster")
    VisReportServer("slide", pgx, output_format="slide")
    AiReportServer("aireport", pgx, save_pgx = save_pgx)
    InfographicServer("infographic", pgx)
    
  }) ## end moduleServer
}
