##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## ---------------------------------------------------------------------
## Shared helpers for the AI Studio generation modules (AiReportServer,
## InfographicServer). Both gate paid generation on the same ownership
## predicate and attribute cost telemetry to the same user; keeping that
## logic here is the single source of truth so the two can't drift.
## ---------------------------------------------------------------------

#' Resolve the current user's email for cost telemetry attribution.
#'
#' @param user_email A reactive/closure over auth$email (read lazily so it
#'   is populated post-login), or a plain value. Falls back to NA.
#' @return Character scalar email, or NA_character_ when unavailable.
studio_user_email <- function(user_email) {
  em <- if (is.function(user_email)) user_email() else user_email
  if (is.null(em) || !nzchar(em)) NA_character_ else em
}

#' Pre-flight authorization gate for paid AI generation.
#'
#' Returns TRUE when the user could persist the result (owner or admin);
#' otherwise shows a styled notification and returns FALSE so callers can
#' fail fast before spending on a generation doomed to no-op on save.
#'
#' @param can_save_pgx Predicate closure over auth/source dir, or NULL (no gate).
#' @param pgx ReactiveValues PGX object.
#' @param session Shiny session used for the notification.
#' @param content Human label for the content type, e.g. "AI reports".
#' @return TRUE if generation is allowed, FALSE otherwise.
studio_can_generate <- function(can_save_pgx, pgx, session,
                                content = "AI content") {
  if (is.null(can_save_pgx)) return(TRUE)
  if (isTRUE(can_save_pgx(pgx))) return(TRUE)
  shiny::showNotification(
    paste0("You can only generate ", content, " for datasets you own."),
    type = "warning", session = session)
  FALSE
}


#' The application server-side logic
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @export
#' Note: pgx needs to be reactiveValues
#'
#'
StudioServer <- function(id, pgx, save_pgx = NULL, can_save_pgx = NULL,
                         user_email = NULL) {
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
    
    observeEvent(input$show_custom, {
      bslib::nav_select("studiopanel","Custom")
      bslib::nav_select("settings","Custom")
    })
    
    ## ----------------- output -------------------------------
    VisReportServer("poster", pgx, output_format="poster")
    VisReportServer("slide", pgx, output_format="slide")
    AiReportServer("aireport", pgx, save_pgx = save_pgx,
        can_save_pgx = can_save_pgx, user_email = user_email)
    InfographicServer("infographic", pgx, save_pgx = save_pgx,
        can_save_pgx = can_save_pgx, user_email = user_email)
    
  }) ## end moduleServer
}
