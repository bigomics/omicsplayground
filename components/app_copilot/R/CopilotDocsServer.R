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
#
# Uploads are gated on extension and size. The agent only ingests
# .pdf/.txt/.md and rejects anything whose extracted text exceeds its own
# limit; a file the UI accepts but the agent drops appears in this list
# while being invisible to the assistant, which is how a user ends up
# believing the chatbot "cannot read" a file it was never given.

# character: extensions the agent can ingest. Must stay in step with
# omicsagentovi's `.ovi_bootstrap_docs_from_dir()` glob and
# `ovi_upload_document()`, which reject everything else.
.COPILOT_DOC_EXTS <- c("pdf", "txt", "md")

# character: the same set as a filename pattern, used to list the directory
# so the panel can never show a file the agent will not pick up.
.COPILOT_DOC_PATTERN <- "\\.(pdf|txt|md)$"

# numeric: per-file upload ceiling in bytes, mirroring the agent's 20,000,000
# extracted-character limit. Bytes are an approximation of characters: for
# text and Markdown the two are within a rounding error, while a PDF extracts
# to far less text than it occupies on disk, so this is conservative for
# PDFs — a large scanned PDF may be refused here even though the agent could
# have parsed it. Refusing predictably beats accepting and silently dropping.
.COPILOT_DOC_MAX_BYTES <- 20000000

#' Classify uploaded files against the agent's ingest rules.
#'
#' Pure helper so the gate is testable without a Shiny session.
#'
#' @param names character vector of uploaded file names.
#' @param sizes numeric vector of file sizes in bytes.
#' @return list(ok = logical, reason = character). `reason` is `NA` for
#'   accepted files and a user-facing explanation otherwise.
#' @noRd
.copilot_filter_uploads <- function(names, sizes) {
  n <- length(names)
  ok     <- rep(TRUE, n)
  reason <- rep(NA_character_, n)
  if (n == 0L) return(list(ok = ok, reason = reason))

  sizes <- suppressWarnings(as.numeric(sizes))
  ext   <- tolower(tools::file_ext(names))

  bad_ext <- !ext %in% .COPILOT_DOC_EXTS
  if (any(bad_ext)) {
    ok[bad_ext] <- FALSE
    label <- ifelse(nzchar(ext[bad_ext]), paste0(".", ext[bad_ext]),
                    "no extension")
    reason[bad_ext] <- sprintf(
      "%s is not supported - upload %s",
      label, paste(paste0(".", .COPILOT_DOC_EXTS), collapse = ", ")
    )
  }

  too_big <- ok & !is.na(sizes) & sizes > .COPILOT_DOC_MAX_BYTES
  if (any(too_big)) {
    ok[too_big] <- FALSE
    reason[too_big] <- sprintf(
      "too large (%.1f MB) - the limit is %.0f MB",
      sizes[too_big] / 1e6, .COPILOT_DOC_MAX_BYTES / 1e6
    )
  }

  empty <- ok & !is.na(sizes) & sizes == 0
  if (any(empty)) {
    ok[empty] <- FALSE
    reason[empty] <- "file is empty"
  }

  list(ok = ok, reason = reason)
}

#' Copilot docs panel UI.
#' @param id Shiny module id.
#' @export
CopilotDocsUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(
    class = "copilot-context-panel",
    style = "padding: 12px 0px;",
    div("Obi knows your current dataset but here you can upload documents as extra context. You can augment or ask Obi to compare findings in the documents."),
    br(),
    br(),    
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
      # Filter on the same pattern the agent globs. Anything else in the
      # directory (a stray .json, a file staged by an older build) would
      # otherwise be listed here as if it were available to the assistant.
      files <- list.files(docs_dir, full.names = TRUE,
                          pattern = .COPILOT_DOC_PATTERN, ignore.case = TRUE)
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
      up <- input$upload

      # Server-side re-check. fileInput's `accept` is a browser hint only —
      # drag-and-drop and several browsers ignore it — so it cannot be the
      # only gate.
      verdict <- .copilot_filter_uploads(up$name, up$size)

      for (i in which(verdict$ok)) {
        file.copy(
          up$datapath[[i]],
          file.path(docs_dir, up$name[[i]]),
          overwrite = TRUE
        )
      }

      rejected <- which(!verdict$ok)
      if (length(rejected)) {
        shiny::showNotification(
          shiny::tags$div(
            shiny::tags$b(sprintf(
              "%d file%s not uploaded",
              length(rejected), if (length(rejected) == 1L) "" else "s"
            )),
            shiny::tags$ul(
              style = "margin: 4px 0 0 0; padding-left: 18px;",
              lapply(rejected, function(i) {
                shiny::tags$li(sprintf("%s - %s", up$name[[i]],
                                       verdict$reason[[i]]))
              })
            )
          ),
          type = "error",
          duration = 12
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
