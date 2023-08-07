##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


shinyalert_storage_full <- function() {
  msg <- paste(
    "You have reached your dataset limit. Please delete some datasets,",
    "<a href='https://events.bigomics.ch/upgrade' target='_blank'>",
    "<b><u>upgrade</u></b></a> your account, or contact us if you want to enter our",
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
