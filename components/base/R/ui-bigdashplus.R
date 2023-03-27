##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## Extra code for bigDash UI
##
##

## This calls the select-bigtab JS in app/R/www/temp.js
bigdash.selectTab <- function(session, selected) {
  shiny:::validate_session_object(session)
  msg <- shiny:::dropNulls(list(value = selected))
  ##session$sendInputMessage(inputId, message)
  session$sendCustomMessage("bigdash-select-tab", msg)
}

bigdash.showTabsGoToDataView <- function(session) {
  session$sendCustomMessage("show-tabs", list())  ## in app/R/www/temp.js
}

## ------------------- sideBar ---------------------------
bigdash.openSidebar <- function() {
  shinyjs::runjs("sidebarOpen()")  ## in app/R/www/temp.js
}

bigdash.closeSidebar <- function() {
  shinyjs::runjs("sidebarClose()")  ## in app/R/www/temp.js
}

bigdash.toggleSidebar <- function(state) {
  if(state) bigdash.openSidebar()
  if(!state) bigdash.closeSidebar()
}

## --------------------menuItem --------------------------
bigdash.showMenuItem <- function(session, item) {
  shiny:::validate_session_object(session)
  msg <- shiny:::dropNulls(list(value = item))
  session$sendCustomMessage("bigdash-show-menuitem", msg)  ## in app/R/www/temp.js
}

bigdash.hideMenuItem <- function(session, item) {
  shiny:::validate_session_object(session)
  msg <- shiny:::dropNulls(list(value = item))
  session$sendCustomMessage("bigdash-hide-menuitem", msg)  ## in app/R/www/temp.js
}

bigdash.toggleMenuItem <- function(session, item, state) {
  if(state) bigdash.showMenuItem(session, item)
  if(!state) bigdash.hideMenuItem(session, item)
}

## --------------------BigTab --------------------------
bigdash.showTab <- function(session, tab) {
  shiny:::validate_session_object(session)
  msg <- shiny:::dropNulls(list(value = tab))
  session$sendCustomMessage("bigdash-show-tab", msg)  ## in app/R/www/temp.js
  bigdash.showMenuItem(session, tab)
}

bigdash.hideTab <- function(session, tab) {
  shiny:::validate_session_object(session)
  msg <- shiny:::dropNulls(list(value = tab))
  session$sendCustomMessage("bigdash-hide-tab", msg)  ## in app/R/www/temp.js
  bigdash.hideMenuItem(session, tab)
}

bigdash.toggleTab <- function(session, tab, state) {
  message("[bigdash.toggleTab] tab = ",tab)
  message("[bigdash.toggleTab] state = ",state)
  if(state) bigdash.showTab(session, tab)
  if(!state) bigdash.hideTab(session, tab)
}
