##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

UploadInputs <- function(id) {
  ns <- shiny::NS(id)  ## namespace
  bigdash::tabSettings(
    shiny::actionLink(ns("module_info"), "Tutorial", icon = shiny::icon("youtube"))
  )
}

UploadUI <- function(id) {
  ns <- shiny::NS(id)  ## namespace
  shiny::tagList(
    UploadModuleUI(ns("upload_panel"))
  )
}
