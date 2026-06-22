# CopilotDocsServer.R — Document upload/list panel for the Copilot board.
#
# Pure-filesystem module: stores uploads in `docs_dir`. Each row carries a
# tickbox that stages the document body into the next user turn (via the
# run controller -> .copilot_stage_user_docs_context -> agent_inject_text_block
# path). Consumed rows are greyed out and disabled until reset_consumed() is
# called by the lifecycle controllers on new chat / tier change / restore.
#
# Delete is a per-row trash icon that writes the file name into a single
# shared input (input$delete_target) so we don't need one observer per row.

#' Copilot docs panel UI.
#' @param id Shiny module id.
#' @export
CopilotDocsUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(
    class = "copilot-context-panel",
    style = "padding: 12px 0px;",
    shiny::fileInput(
      ns("upload"),
      label = NULL,
      accept = c(".txt", ".pdf", ".md"),
      buttonLabel = "Upload...",
      multiple = TRUE,
      width = "100%"
    ),
    shiny::uiOutput(ns("docs_ui"))
  )
}

#' Copilot docs panel server.
#'
#' @param id Shiny module id.
#' @param docs_dir character(1) path to the docs directory. Must exist
#'   before the module is constructed.
#' @return list(selected_docs, mark_consumed, reset_consumed, doc_count,
#'   on_delete).
#' @export
CopilotDocsServer <- function(id, docs_dir) {
  shiny::moduleServer(id, function(input, output, session) {

    .refresh   <- shiny::reactiveVal(0L)
    .on_delete <- shiny::reactiveVal(NULL)
    consumed   <- shiny::reactiveVal(character(0))

    .doc_files <- shiny::reactive({
      .refresh()
      files <- list.files(docs_dir, full.names = TRUE)
      if (!length(files)) {
        return(data.frame(Name = character(0), Size = character(0),
                          stringsAsFactors = FALSE))
      }
      info <- file.info(files)
      data.frame(
        Name = basename(files),
        Size = paste0(round(info$size / 1024, 1), " KB"),
        stringsAsFactors = FALSE
      )
    })

    .input_id <- function(i) paste0("doc_", i)

    shiny::observeEvent(input$upload, {
      shiny::req(input$upload)
      for (i in seq_len(nrow(input$upload))) {
        file.copy(
          input$upload$datapath[[i]],
          file.path(docs_dir, input$upload$name[[i]]),
          overwrite = TRUE
        )
      }
      .refresh(.refresh() + 1L)
    })

    shiny::observeEvent(input$delete_target, {
      nm <- input$delete_target
      shiny::req(nm)
      fpath <- file.path(docs_dir, nm)
      if (file.exists(fpath)) file.remove(fpath)
      .on_delete(fpath)
      consumed(setdiff(shiny::isolate(consumed()), nm))
      .refresh(shiny::isolate(.refresh()) + 1L)
    }, ignoreInit = TRUE)

    selected_docs <- shiny::reactive({
      df <- .doc_files()
      if (!nrow(df)) return(character(0))
      used <- intersect(consumed(), df$Name)
      keep <- vapply(seq_len(nrow(df)), function(i) {
        nm <- df$Name[[i]]
        if (nm %in% used) return(FALSE)
        val <- input[[.input_id(i)]]
        isTRUE(val)
      }, logical(1))
      df$Name[keep]
    })

    mark_consumed <- function(docs) {
      docs <- tryCatch(as.character(docs), error = function(e) character(0))
      docs <- docs[!is.na(docs) & nzchar(docs)]
      if (!length(docs)) return(invisible(NULL))
      consumed(unique(c(shiny::isolate(consumed()), docs)))
      invisible(NULL)
    }

    reset_consumed <- function() {
      consumed(character(0))
      invisible(NULL)
    }

    output$docs_ui <- shiny::renderUI({
      df <- .doc_files()
      if (!nrow(df)) {
        return(shiny::tags$div(
          class = "text-muted small",
          "No documents uploaded."
        ))
      }
      used <- intersect(consumed(), df$Name)
      shiny::tags$div(
        shiny::tags$div(
          class = "small text-muted",
          style = "margin-bottom: 4px;",
          "Uploaded documents"
        ),
        shiny::tags$ul(
          style = paste(
            "list-style: none;",
            "padding-left: 0;",
            "margin: 0;"
          ),
          lapply(seq_len(nrow(df)), function(i) {
            nm <- df$Name[[i]]
            sz <- df$Size[[i]]
            id <- .input_id(i)
            already_used <- nm %in% used
            current <- input[[id]]
            checked <- if (already_used) {
              TRUE
            } else if (is.null(current)) {
              FALSE
            } else {
              isTRUE(current)
            }
            shiny::tags$li(
              style = paste(
                "display: flex; align-items: center;",
                "padding: 2px 0;",
                if (already_used) "opacity: 0.55;" else ""
              ),
              shiny::tags$label(
                class = "checkbox-inline",
                style = paste(
                  "display: flex; align-items: center;",
                  "font-weight: normal; flex: 1; margin: 0;"
                ),
                shiny::tags$input(
                  id = session$ns(id),
                  type = "checkbox",
                  checked = if (checked) "checked" else NULL,
                  disabled = if (already_used) "disabled" else NULL
                ),
                shiny::tags$span(
                  style = "margin-left: 6px; word-break: break-all;",
                  nm
                ),
                shiny::tags$span(
                  class = "small text-muted",
                  style = "margin-left: 8px;",
                  sz
                )
              ),
              shiny::tags$button(
                type = "button",
                class = "btn btn-link btn-sm text-danger",
                style = "padding: 0 4px; margin-left: 4px;",
                title = "Delete",
                onclick = sprintf(
                  "Shiny.setInputValue('%s', '%s', {priority: 'event'})",
                  session$ns("delete_target"), nm
                ),
                shiny::tags$i(class = "fa fa-trash")
              )
            )
          })
        )
      )
    })

    list(
      selected_docs  = selected_docs,
      mark_consumed  = mark_consumed,
      reset_consumed = reset_consumed,
      doc_count      = shiny::reactive(nrow(.doc_files())),
      on_delete      = .on_delete
    )
  })
}
