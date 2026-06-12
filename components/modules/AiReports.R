##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# PGX AI report access helpers
# =============================================================================

.ai_report_ai_slot <- function(pgx) {
  if (is.null(pgx) || !is.list(pgx)) return(NULL)
  ai <- pgx$ai
  if (is.null(ai) || !is.list(ai)) return(NULL)
  ai
}

.ai_report_valid_entry <- function(x) {
  is.list(x) &&
    is.character(x$report) &&
    length(x$report) > 0L &&
    !is.na(x$report[[1L]]) &&
    nzchar(x$report[[1L]])
}

ai_report_slots <- function(pgx) {
  ai <- .ai_report_ai_slot(pgx)
  if (is.null(ai)) return(character(0))

  slots <- setdiff(names(ai), "meta")
  slots[vapply(slots, function(slot) {
    .ai_report_valid_entry(ai[[slot]])
  }, logical(1))]
}

ai_report_get <- function(pgx, slot) {
  if (is.null(slot) || length(slot) != 1L) return(NULL)
  slot <- tryCatch(as.character(slot)[[1L]], error = function(e) NA_character_)
  if (is.na(slot) || !nzchar(slot)) {
    return(NULL)
  }

  ai <- .ai_report_ai_slot(pgx)
  if (is.null(ai)) return(NULL)

  entry <- ai[[slot]]
  if (!.ai_report_valid_entry(entry)) return(NULL)

  list(
    slot = slot,
    report = as.character(entry$report)[[1L]],
    prompt = if (is.character(entry$prompt) && length(entry$prompt) > 0L) {
      as.character(entry$prompt)[[1L]]
    } else {
      NULL
    },
    meta = if (is.null(ai$meta)) NULL else ai$meta
  )
}

ai_report_has <- function(pgx, select = NULL) {
  slots <- ai_report_slots(pgx)
  if (is.null(select)) return(length(slots) > 0L)
  select <- tryCatch(as.character(select), error = function(e) character(0))
  select <- select[!is.na(select) & nzchar(select)]
  if (!length(select)) return(FALSE)
  all(vapply(select, function(slot) {
    slot %in% slots ||
      (identical(slot, "drugs") && any(startsWith(slots, "drugs_")))
  }, logical(1)))
}

ai_report_needs_generation <- function(pgx) {
  length(ai_report_modules_for_pgx(pgx)) > 0L && !ai_report_has(pgx)
}

ai_report_drug_slots <- function(pgx) {
  slots <- ai_report_slots(pgx)
  slots[startsWith(slots, "drugs_")]
}

ai_report_drug_label <- function(pgx, slot) {
  slot <- tryCatch(as.character(slot)[[1L]], error = function(e) NA_character_)
  if (is.na(slot) || !nzchar(slot)) return("")

  label <- gsub("_", " ", sub("^drugs_", "", slot))
  dbs <- tryCatch(names(pgx$drugs), error = function(e) character(0))
  if (length(dbs)) {
    safe <- paste0("drugs_", gsub("[^A-Za-z0-9]+", "_", dbs))
    match_idx <- match(slot, safe)
    if (!is.na(match_idx)) label <- dbs[[match_idx]]
  }

  if (grepl("^L1000([_/ -]|$)", label, ignore.case = TRUE)) {
    if (grepl("activ", label, ignore.case = TRUE)) return("L1000 Activity")
    if (grepl("gene", label, ignore.case = TRUE)) return("L1000 Gene")
  }
  label
}

ai_report_modules_for_pgx <- function(pgx) {
  if (is.null(pgx) || !is.list(pgx)) return(character(0))

  modules <- c(character(0),
    if (!is.null(pgx$wgcna)) "wgcna",
    if (!is.null(pgx$wgcna_mox)) "wgcna_mox",
    if (!is.null(pgx$mofa)) "mofa",
    if (!is.null(pgx$drugs) && length(pgx$drugs) > 0L) "drugs",
    if (!is.null(pgx$gx.meta)) "de",
    if (!is.null(pgx$gset.meta)) "pathways"
  )
  if (length(modules) > 0L) modules <- c(modules, "combined")
  unique(modules)
}

ai_report_get_module <- function(pgx, module) {
  if (is.null(module) || length(module) != 1L) return(NULL)
  module <- tryCatch(
    as.character(module)[[1L]],
    error = function(e) NA_character_
  )
  if (is.na(module) || !nzchar(module)) return(NULL)

  module <- switch(module,
    summary = "combined",
    enrichment = "pathways",
    module
  )
  ai_report_get(pgx, module)
}

ai_report_dataset_token <- function(pgx) {
  if (is.null(pgx) || !is.list(pgx)) return("")
  x <- tryCatch(pgx$X, error = function(e) NULL)
  xdim <- if (is.null(dim(x))) "" else paste(dim(x), collapse = "x")
  name <- if (is.null(pgx$name)) "" else as.character(pgx$name)[[1L]]
  date <- if (is.null(pgx$date)) "" else as.character(pgx$date)[[1L]]
  paste(name, xdim, date, sep = "|")
}

ai_report_generate <- function(pgx,
                               llm_model,
                               force = FALSE,
                               select = NULL,
                               img_model = NULL,
                               report_type = "normal",
                               on_error = "warn") {
  if (is.null(llm_model) || !nzchar(llm_model)) return(pgx)
  if (is.null(select)) select <- ai_report_modules_for_pgx(pgx)
  select <- tryCatch(as.character(select), error = function(e) character(0))
  select <- select[!is.na(select) & nzchar(select)]
  if (!length(select)) return(pgx)

  playbase::pgx.update_reports(
    pgx,
    ai = list(
      llm_model = llm_model,
      img_model = img_model,
      select = select,
      report_type = report_type,
      force = isTRUE(force),
      on_error = on_error
    )
  )
}

ai_report_update_text <- function(pgx, reports) {
  if (is.null(pgx) || !is.list(pgx)) return(pgx)
  ai <- .ai_report_ai_slot(pgx)
  if (is.null(ai)) return(pgx)

  reports <- reports[!is.na(names(reports)) & nzchar(names(reports))]
  if (!length(reports)) return(pgx)

  for (slot in names(reports)) {
    if (!is.list(ai[[slot]])) next
    value <- tryCatch(as.character(reports[[slot]])[[1L]],
      error = function(e) NA_character_)
    if (is.na(value)) next
    ai[[slot]]$report <- value
  }

  pgx$ai <- ai
  pgx
}

ai_report_merge_into_reactive <- function(pgx_rv, ai) {
  if (is.null(pgx_rv)) return(invisible(FALSE))

  ai_slot <- if (is.list(ai) && !is.null(ai$ai) && is.list(ai$ai)) ai$ai else ai
  probe <- list(ai = ai_slot)
  slots <- ai_report_slots(probe)
  if (is.null(ai_slot) || !is.list(ai_slot) || !length(slots)) {
    return(invisible(FALSE))
  }

  current <- pgx_rv$ai
  if (is.null(current) || !is.list(current)) current <- list()
  for (slot in slots) {
    current[[slot]] <- ai_slot[[slot]]
  }
  if (!is.null(ai_slot$meta)) current$meta <- ai_slot$meta
  pgx_rv$ai <- current
  invisible(TRUE)
}

ai_report_copy_into_reactive <- function(pgx_rv, ai) {
  if (is.null(pgx_rv)) return(invisible(FALSE))

  ai_slot <- if (is.list(ai) && !is.null(ai$ai) && is.list(ai$ai)) ai$ai else ai
  probe <- list(ai = ai_slot)
  if (is.null(ai_slot) || !is.list(ai_slot) || !ai_report_has(probe)) {
    return(invisible(FALSE))
  }

  pgx_rv$ai <- ai_slot
  invisible(TRUE)
}

#' Get one stored AI infographic from pgx$ai
#'
#' Returns an infographic entry only when it contains stored bytes, an existing
#' file path, or a stored error status.
#'
#' @param pgx PGX object as a list.
#' @param slot AI report slot name.
#' @return Infographic entry list, or \code{NULL} when unavailable.
ai_infographic_get <- function(pgx, slot) {
  if (is.null(slot) || length(slot) != 1L) return(NULL)
  slot <- tryCatch(as.character(slot)[[1L]], error = function(e) NA_character_)
  if (is.na(slot) || !nzchar(slot)) return(NULL)

  ai <- .ai_report_ai_slot(pgx)
  if (is.null(ai) || !is.list(ai[[slot]])) return(NULL)

  img <- ai[[slot]]$infographic
  if (!is.list(img)) return(NULL)
  has_bytes <- is.raw(img$bytes) && length(img$bytes) > 0L
  has_path <- is.character(img$path) && length(img$path) > 0L &&
    nzchar(img$path[[1L]]) && file.exists(img$path[[1L]])
  if (!has_bytes && !has_path && !identical(img$status, "error")) return(NULL)
  img
}

#' List AI report slots with stored infographic state
#'
#' @param pgx PGX object as a list.
#' @return Character vector of AI report slot names with infographic entries.
ai_infographic_slots <- function(pgx) {
  ai <- .ai_report_ai_slot(pgx)
  if (is.null(ai)) return(character(0))

  slots <- setdiff(names(ai), "meta")
  slots[vapply(slots, function(slot) {
    !is.null(ai_infographic_get(pgx, slot))
  }, logical(1))]
}

#' Convert provider errors into user-facing infographic messages
#'
#' @param error Error condition or character error message.
#' @return Sanitized message suitable for storage and UI display.
ai_infographic_friendly_error <- function(error) {
  msg <- if (inherits(error, "condition")) {
    conditionMessage(error)
  } else if (is.null(error)) {
    ""
  } else {
    paste(as.character(error), collapse = " ")
  }
  msg <- trimws(msg)

  if (!nzchar(msg)) {
    return("Infographic generation failed. Please try again.")
  }
  if (grepl("503|Service Unavailable|high demand|overload|saturat",
      msg, ignore.case = TRUE)) {
    return(paste0(
      "Image generation server is temporarily overloaded. ",
      "Please try again in a few minutes."
    ))
  }
  if (grepl("429|Too Many Requests|rate.?limit", msg, ignore.case = TRUE)) {
    return("Image generation rate limit reached. Please wait a moment and try again.")
  }
  if (grepl("timeout|timed.?out", msg, ignore.case = TRUE)) {
    return(paste0(
      "Image generation timed out. The service may be busy; ",
      "please try again."
    ))
  }
  if (grepl("No image data|Empty response|no response", msg, ignore.case = TRUE)) {
    return(paste0(
      "Image generation server seems saturated and returned no image. ",
      "Please try again."
    ))
  }
  if (grepl("_API_KEY|API.?key|401|Unauthorized", msg, ignore.case = TRUE)) {
    return(paste0(
      "Image generation is not configured correctly. ",
      "Please contact your administrator."
    ))
  }

  "Image generation failed. Please try again."
}

#' Extract persisted image payload fields from an omicsai image result
#'
#' @param result \code{omicsai_image_result} list, or \code{NULL} for errors.
#' @return List containing bytes, content type, path, prompt, model, and metadata.
ai_infographic_payload <- function(result) {
  path <- NULL
  prompt <- NULL
  metadata <- NULL
  if (is.list(result)) {
    path <- if (is.character(result$path) && length(result$path) > 0L) {
      result$path[[1L]]
    } else {
      NULL
    }
    prompt <- result$prompt
    metadata <- result$metadata
  }

  bytes <- raw(0)
  if (!is.null(path) && file.exists(path)) {
    bytes <- readBin(path, what = raw(), n = file.info(path)$size)
  }

  ext <- tolower(tools::file_ext(if (is.null(path)) "image.png" else path))
  list(
    bytes = bytes,
    content_type = switch(ext,
      jpg = "image/jpeg",
      jpeg = "image/jpeg",
      "image/png"),
    path = path,
    prompt = prompt,
    model = if (is.null(metadata$model)) NULL else metadata$model,
    metadata = metadata
  )
}

#' Store one AI infographic under pgx$ai
#'
#' Writes generated image bytes and metadata, or a sanitized error entry, into
#' \code{pgx$ai[[slot]]$infographic}.
#'
#' @param pgx PGX object as a list.
#' @param slot AI report slot name.
#' @param result \code{omicsai_image_result} list, or \code{NULL} for errors.
#' @param status Infographic status, usually \code{"done"} or \code{"error"}.
#' @param error Optional provider error for failed generation.
#' @param style omicsai image style ID.
#' @param n_blocks Deprecated layout metadata retained for compatibility.
#' @return Updated PGX object.
ai_infographic_set <- function(pgx, slot, result,
                               status = "done", error = NULL,
                               style = NULL, n_blocks = NULL) {
  if (is.null(pgx) || !is.list(pgx)) return(pgx)
  slot <- tryCatch(as.character(slot)[[1L]], error = function(e) NA_character_)
  if (is.na(slot) || !nzchar(slot)) return(pgx)

  if (is.null(pgx$ai) || !is.list(pgx$ai)) pgx$ai <- list()
  if (is.null(pgx$ai[[slot]]) || !is.list(pgx$ai[[slot]])) {
    pgx$ai[[slot]] <- list()
  }
  # Normalize generated image data before writing the PGX AI schema entry.
  payload <- ai_infographic_payload(result)

  pgx$ai[[slot]]$infographic <- list(
    status = status,
    error = if (identical(status, "error")) {
      ai_infographic_friendly_error(error)
    } else {
      error
    },
    bytes = payload$bytes,
    content_type = payload$content_type,
    path = payload$path,
    prompt = payload$prompt,
    model = payload$model,
    style = style,
    n_blocks = n_blocks,
    date = as.character(Sys.time()),
    metadata = payload$metadata
  )
  pgx
}

#' Convert a stored AI infographic entry into a Shiny renderImage value
#'
#' @param img Infographic entry returned by \code{ai_infographic_get()}.
#' @param tmpdir Directory used to materialize stored raw bytes.
#' @param name Base filename for materialized image bytes.
#' @return List suitable for \code{shiny::renderImage()}, or \code{NULL}.
ai_infographic_render_value <- function(img, tmpdir, name = "infographic") {
  if (is.null(img) || !is.list(img)) return(NULL)
  if (identical(img$status, "error")) {
    # Keep legacy/raw stored errors friendly before they reach the UI.
    msg <- ai_infographic_friendly_error(img$error)
    shiny::validate(shiny::need(FALSE, msg))
  }

  content_type <- if (is.null(img$content_type)) "image/png" else img$content_type
  ext <- if (identical(content_type, "image/jpeg")) ".jpg" else ".png"
  src <- NULL
  if (is.raw(img$bytes) && length(img$bytes) > 0L) {
    dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
    src <- file.path(tmpdir, paste0(name, ext))
    writeBin(img$bytes, src)
  } else if (is.character(img$path) && length(img$path) > 0L &&
      file.exists(img$path[[1L]])) {
    src <- img$path[[1L]]
  }

  if (is.null(src)) return(NULL)
  list(src = src, height = "auto", width = "100%", contentType = content_type)
}
