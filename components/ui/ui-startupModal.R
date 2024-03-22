##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ui.startupModal <- function(id, messages, title) {
  if (length(messages) == 0) {
    return(NULL)
  }

  carousel_items <- list()
  for (i in 1:length(messages)) {
    tag1 <- bsutils::carouselItem(
      div(
        style = "height: 360px;",
        class = "d-flex align-items-center justify-content-center",
        HTML(paste0("<div>", messages[[i]], "</div>"))
      ),
      class = "p-4"
    )
    carousel_items[[i]] <- tag1
  }

  modal <- shiny::modalDialog(
    size = "l",
    title = NULL,
    footer = NULL,
    bsutils::modalHeader(
      bsutils::modalTitle(title),
      style = "background-color: #f0f9fd;"),
    do.call(
      function(...) {
        bsutils::carousel(
          ...,
          id = "opg-welcome-carousel",
          indicators = TRUE,
          controls = TRUE
        )
      },
      carousel_items
    ),
    easyClose = TRUE,
    tags$style(".modal-dialog {width: 720px;}"),
    tags$style(".modal-content {background-color: #f0f9fd;}"),    
    tags$style(".modal-header {padding: 0px;}")    
  )
  ##return(div(id=id, modal))
  return(modal)
}


ui.inviteModal <- function(id) {

  body <- tagList(
    HTML("<center><h3>and earn some BigOmics swag!</h3><p><p>Invite your friends to Omics Playground and earn a cool Bigomics 'Friendly Monster' T-Shirt (minimum 3 referrals) or one of our awesome sustainable Dopper water bottles (minimum 10 referrals)."),
    HTML("<br><br><p>Enter your friend's email:"),
    div( shiny::textInput(
        inputId = paste0(id,"_email"),
        label = "", placeholder = "Email address..."),
        style = "margin-top: -30px;"),
    br(),
    shiny::actionButton(paste0(id,"_button"),"Invite!", class="btn btn-primary" )
  )          

  modal <- shiny::modalDialog(
    title = NULL,
    bsutils::modalHeader(
      div(class="modal-title", "Invite A Friend"),
      style = "background-color: #f0f9fd;"),
    body,
    footer = NULL,
    size = "l",
    easyClose = TRUE,
    tags$style(".modal-dialog {width: 720px;}"),
    tags$style(".modal-content {background-color: #f0f9fd;}"),    
    tags$style(".modal-header {padding: 0px;}")    
  )
  ##return(div(id=id, modal))
  return(shiny::showModal(modal))
}
