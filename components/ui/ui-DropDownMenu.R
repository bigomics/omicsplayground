##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

DropdownMenu <- function(..., size = "default", status = "default", icon = NULL, circle = TRUE, border = "default", width = "") {
  id <- bigdash:::make_id()
  bslib::popover(
    tags$a(
      class = paste0(
        "btn btn-", status, " action-button",
        if (circle) {
          " btn-circle"
        } else {
          ""
        },
        if (size == "default") {
          " "
        } else {
          paste0("-", size, " ")
        },
        if (size != "default") paste0("btn-", size),
        if (border != "default") paste0(" border border-", border)
      ),
      list(icon), `data-bs-toggle` = "popover", id = id
    ),
    ...,
    placement = "bottom",
    options = list(
      animation = FALSE,
      template = paste0('<div class="popover" role="tooltip"><div class="popover-arrow"></div><h3 class="popover-header"></h3><div class="popover-body style=width:,', width, ';"></div></div>')
    ),
    onClick = HTML(paste0("$(document).ready(function() { makePopoverDismissible('", id, "'); addActionOnPopoverChange('", id, "'); });"))
  )
}
