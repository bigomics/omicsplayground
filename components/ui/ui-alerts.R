##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


shinyalert_storage_full <- function(numpgx = NULL, maxpgx = NULL) {
  msg <- paste("You have reached your dataset limit.")
  if (!is.null(numpgx) && !is.null(maxpgx)) {
    msg <- paste(
      "You currently have", numpgx, "datasets.",
      "You have reached your dataset limit of", maxpgx, "."
    )
  }
  msg <- paste(
    msg, "Please <a href='https://events.bigomics.ch/upgrade' target='_blank'>",
    "<b><u>upgrade</u></b></a> or get free extra datasets by entering our",
    "<a href='https://bigomics.ch/contact-us' target='_blank'><b><u>Ambassador",
    "or Promoter</u></b></a> program."
  )

  shinyalert::shinyalert(
    title = "Your storage is full!",
    text = HTML(msg),
    html = TRUE,
    type = "warning"
  )
}
