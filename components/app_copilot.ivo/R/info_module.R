copilot_info_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("dataset_info"))
}


copilot_info_module <- function(id, pgx) {
  shiny::moduleServer(id, function(input, output, session) {

    # ---- Dataset context card output ----
    output$dataset_info <- shiny::renderUI({
      shiny::validate(shiny::need(
        !is.null(pgx) && !is.null(pgx$X),
        "No dataset loaded yet."
      ))
      ##pgx <- local_pgx()

      n_samples <- if (!is.null(pgx$X)) ncol(pgx$X) else "?"
      n_genes   <- if (!is.null(pgx$X)) nrow(pgx$X) else "?"
      organism  <- if (!is.null(pgx$organism)) pgx$organism else "unknown"
      dataname  <- if (!is.null(pgx$name)) pgx$name else "unknown"
      description  <- if (!is.null(pgx$description)) pgx$description else "none"            

      contrasts <- "none"
      if (!is.null(pgx$contrasts)) {
        ct <- colnames(pgx$contrasts)
        if (length(ct) > 3L) {
          contrasts <- paste(c(ct[1:3], paste0("+ ", length(ct) - 3L, " more")), collapse = ", ")
        } else {
          contrasts <- paste(ct, collapse = ", ")
        }
      }

      shiny::tags$div(
        style = "font-size: 0.9em; line-height: 1.6em; padding: 1.3em 0.6em;",
        shiny::tags$div(shiny::strong("Dataset: "),  dataname),
        shiny::tags$div(shiny::strong("Description: "),  description),
        shiny::tags$div(shiny::strong("Organism: "),  organism),                
        shiny::tags$div(shiny::strong("Samples: "),   n_samples),
        shiny::tags$div(shiny::strong("Genes: "),     n_genes),
        shiny::tags$div(shiny::strong("Contrasts: "), contrasts)
      )
    })

    
  })
}
