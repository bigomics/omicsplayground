##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

BIGOMICS_CONTACT_US_URL <- "https://bigomics.ch/contact-us/"

contact_us_callback_js <- function() {
  sprintf(
    "function(value) { if (value === true) { window.open('%s', '_blank'); } }",
    BIGOMICS_CONTACT_US_URL
  )
}

## user_level: auth$level from AuthenticationModule ("free", "starter", "premium", ...)
auth_user_level_is_free <- function(user_level) {
  ul <- if (is.null(user_level) || (length(user_level) == 1L && is.na(user_level))) {
    ""
  } else {
    as.character(user_level)
  }
  identical(tolower(trimws(ul)), "free")
}

## Standard Ok (cancel) + Contact us (confirm, opens contact page in new tab)
shinyalert_ok_contact_us <- function(title,
                                     text,
                                     html = FALSE,
                                     type = "warning") {
  shinyalert::shinyalert(
    title = title,
    text = if (html) shiny::HTML(text) else text,
    html = html,
    type = type,
    showCancelButton = TRUE,
    confirmButtonText = "Contact us",
    cancelButtonText = "Ok",
    confirmButtonCol = "#337ab7",
    callbackJS = contact_us_callback_js()
  )
}

## user_level: auth$level from AuthenticationModule ("free", "starter", "premium", ...)
shinyalert_storage_full <- function(numpgx = NULL, maxpgx = NULL, user_level = NULL) {
  is_free <- auth_user_level_is_free(user_level)
  limit_str <- if (!is.null(maxpgx) && length(maxpgx) && !is.na(maxpgx)) {
    as.character(maxpgx)
  } else {
    "X"
  }

  if (is_free) {
    msg <- paste0(
      "You've reached your free dataset limit of ", limit_str,
      ". To upload more, contact our team to upgrade your account or explore our ",
      "<a href='https://events.bigomics.ch/upgrade' target='_blank'>subscription plans here</a>."
    )
    shinyalert::shinyalert(
      title = "Your storage is full!",
      text = shiny::HTML(msg),
      html = TRUE,
      type = "warning",
      confirmButtonText = "Contact us",
      confirmButtonCol = "#337ab7",
      callbackJS = contact_us_callback_js()
    )
  } else {
    msg <- paste(
      "You've reached your limit for concurrent datasets.",
      "Delete one of your existing datasets or contact support to increase your storage."
    )
    shinyalert_ok_contact_us(
      title = "Your storage is full!",
      text = msg,
      type = "warning"
    )
  }
}

## file_kind: "counts" = counts.csv check; "samples" = samples.csv row check
shinyalert_max_samples_reached <- function(max_samples,
                                           user_level = NULL,
                                           file_kind = c("counts", "samples")) {
  file_kind <- match.arg(file_kind)
  file_phrase <- if (file_kind == "counts") "counts file" else "samples file"
  x_str <- if (!is.null(max_samples) && length(max_samples) && !is.na(max_samples)) {
    as.character(max_samples)
  } else {
    "X"
  }
  is_free <- auth_user_level_is_free(user_level)
  if (is_free) {
    msg <- paste0(
      "You have reached the maximum number of samples allowed in the free trial. ",
      "Please upload a new ", file_phrase, " with a maximum of ", x_str,
      " samples or contact us to upgrade your account."
    )
  } else {
    msg <- paste0(
      "You have reached the maximum number of samples allowed for your account. ",
      "Please upload a new ", file_phrase, " with a maximum of ", x_str,
      " samples or contact support."
    )
  }
  shinyalert_ok_contact_us(
    title = "Maximum samples reached",
    text = msg,
    type = "error"
  )
}
