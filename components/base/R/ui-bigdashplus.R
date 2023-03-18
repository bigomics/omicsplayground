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


#--------------------------------------------------------------------------------
# This JS code is responsible for activating and desactivating the class `active`,
# defined on `scss/components/_dropdownmenu.scss` when clicking a dropdown button.
# It also has logic implemented so that when clicking a different dropdown item,
# the active one gets desactivated.
jsCode <- function(id){
  paste0(
    "$('#", id, "').on('click', function(){
    if($('#", id, "').hasClass('active')){
         $('#", id,"').toggleClass('active');
    return 0;
      }
      $('.dropdown-button.active').toggleClass('active');
      if(!$('#", id, "').hasClass('active')){
         $('#", id,"').toggleClass('active');
      }

})"
  )
}

DropdownMenu <- function (..., size = "default", status = "default", icon = NULL, width = "250px", margin = "10px")
{
  id <- bigdash:::make_id()
  tags$div(tags$a(class = paste0(
    "btn btn-", status, " action-button dropdown-button btn-circle", if(size == "default"){" "}else{paste0("-", size, " ")},
    if (size != "default") paste0("btn-", size)
  ),
  id = id, type = "button", `data-bs-toggle` = "dropdown",
  `aria-expanded` = "false", list(icon)),
  tags$div(class = "dropdown-menu",
           `aria-labelledby` = id,
           tryCatch({
             lapply(list(...), function(x){
               htmltools::tags$li(
                 x,
                 style = paste0(
                   "width: ", htmltools::validateCssUnit(width),
                   "; margin-left: ", htmltools::validateCssUnit(margin),                   
                   "; margin-right: ", htmltools::validateCssUnit(margin), ";")
               )
             })
           }, error = function(w){})
  ),
  tags$script(HTML(jsCode(id)))
  )
}

# same thing as DropdownMenu but without circle icon;
# currently used in the board.loading table
actionMenu <- function (..., size = "default", status = "default", icon = NULL, margin = "10px")
{
  id <- bigdash:::make_id()
  tags$div(tags$a(class = paste0(
    "btn btn-outline-", status, " action-button dropdown-button", if(size == "default"){" "}else{paste0("-", size, " ")},
    if (size != "default") paste0("btn-", size)
  ),
  style = 'border: none;',
  id = id, type = "button", `data-bs-toggle` = "dropdown",
  `aria-expanded` = "false", list(icon)),
  tags$div(class = "dropdown-menu",
           `aria-labelledby` = id,
           tryCatch({
             lapply(list(...), function(x){
               htmltools::tags$li(
                 x,
                 style = paste0("margin-left: ",
                                htmltools::validateCssUnit(margin), "; margin-right: ",
                                htmltools::validateCssUnit(margin), ";")
               )
             })
           }, error = function(w){})
  ),
  tags$script(HTML(jsCode(id)))
  )
}
