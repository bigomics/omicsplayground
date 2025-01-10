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
  #
  session$sendCustomMessage("bigdash-select-tab", msg)
}

bigdash.showTabsGoToDataView <- function(session) {
  session$sendCustomMessage("show-tabs", list()) ## in app/R/www/temp.js
}

## ------------------- sideBar ---------------------------
bigdash.openSidebar <- function() {
  shinyjs::runjs("sidebarOpen()") ## in app/R/www/temp.js
}

bigdash.closeSidebar <- function() {
  shinyjs::runjs("sidebarClose()") ## in app/R/www/temp.js
}

bigdash.unloadSidebar <- function() {
  shinyjs::runjs("unloadSidebar()") ## in app/R/www/temp.js
}

bigdash.toggleSidebar <- function(state) {
  if (state) bigdash.openSidebar()
  if (!state) bigdash.closeSidebar()
}

## ------------------- settingsBar ---------------------------
bigdash.openSettings <- function(lock = TRUE) {
  shinyjs::runjs("settingsUnlock()") ## in app/R/www/temp.js
  shinyjs::runjs("settingsOpen()") ## in app/R/www/temp.js
  if (lock) {
    Sys.sleep(0.1)
    shinyjs::runjs("settingsLock()") ## in app/R/www/temp.js
  }
}

bigdash.closeSettings <- function() {
  shinyjs::runjs("settingsUnlock()") ## in app/R/www/temp.js
  shinyjs::runjs("settingsClose()") ## in app/R/www/temp.js
}

## --------------------menuItem --------------------------
bigdash.showMenuItem <- function(session, item) {
  shiny:::validate_session_object(session)
  msg <- shiny:::dropNulls(list(value = item))
  session$sendCustomMessage("bigdash-show-menuitem", msg) ## in app/R/www/temp.js
}

bigdash.hideMenuItem <- function(session, item) {
  shiny:::validate_session_object(session)
  msg <- shiny:::dropNulls(list(value = item))
  session$sendCustomMessage("bigdash-hide-menuitem", msg) ## in app/R/www/temp.js
}

bigdash.toggleMenuItem <- function(session, item, state) {
  if (state) bigdash.showMenuItem(session, item)
  if (!state) bigdash.hideMenuItem(session, item)
}

## --------------------BigTab --------------------------
bigdash.showTab <- function(session, tab) {
  shiny:::validate_session_object(session)
  msg <- shiny:::dropNulls(list(value = tab))
  session$sendCustomMessage("bigdash-show-tab", msg) ## in app/R/www/temp.js
  bigdash.showMenuItem(session, tab)
}

bigdash.hideTab <- function(session, tab) {
  shiny:::validate_session_object(session)
  msg <- shiny:::dropNulls(list(value = tab))
  session$sendCustomMessage("bigdash-hide-tab", msg) ## in app/R/www/temp.js
  bigdash.hideMenuItem(session, tab)
}

bigdash.toggleTab <- function(session, tab, state) {
  if (state) bigdash.showTab(session, tab)
  if (!state) bigdash.hideTab(session, tab)
}

bigdash.removeTab <- function(session, tab) {
  shiny:::validate_session_object(session)
  msg <- shiny:::dropNulls(list(value = tab))
  session$sendCustomMessage("bigdash-remove-tab", msg)
}

bigdash.hideMenuElement <- function(session, name) {
  shiny:::validate_session_object(session)
  msg <- shiny:::dropNulls(list(value = name))
  session$sendCustomMessage("bigdash-hide-menu-element", msg)
}

bigdash.showMenuElement <- function(session, name) {
  shiny:::validate_session_object(session)
  msg <- shiny:::dropNulls(list(value = name))
  session$sendCustomMessage("bigdash-show-menu-element", msg)
}
