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
