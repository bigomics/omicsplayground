##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

gadgetize <- function(moduleUI, moduleSERVER, title="shiny gadget", ...)
{
    ## Creates gadget from a Shiny module. Gadget are browser-based UI
    ## applets for a single task that can be run from the R command
    ## line. It is not to be used inside Shiny programs themselves.
    ##
    
    id = sub(".*file","gadget",tempfile())  ## random ID
    ui = miniUI::miniPage(
        shiny::tags$head(shiny::tags$style(".modal-dialog{width:900px}")),
        shiny::tags$head(shiny::tags$style(".modal-dialog.modal-lg{width:1400px}")),
        shiny::tags$head(shiny::tags$style(".modal-dialog.modal-sm{width:400px}")),        
        ##shinyalert::useShinyalert(),
        miniUI::gadgetTitleBar(title),
        miniUI::miniContentPanel(moduleUI(id))
    )
    server = function(input, output, session) {
        return_obj <- moduleSERVER(id, ...)
        ##return_obj <- moduleSERVER(id)    
        shiny::observeEvent( input$done, {
            shiny::stopApp(return_obj())
        })
    }
    X <- shiny::runGadget(ui, server)
    cat("[gadgetize] names(X)=",names(X),"\n")
    cat("[gadgetize] *** closing gadget ***\n")
    X
}

gadgetize2 <- function(moduleUI, moduleSERVER, title="shiny gadget",
                       size="m", ...)
{
    ##
    ## Creates modalDialog from a Shiny module similar as used inside
    ## Shiny programs.
    ##
    
    
    id = sub(".*file","gadget",tempfile())  ## random ID
    ui = shiny::fluidPage(
        shiny::tags$head(shiny::tags$style(".modal-dialog{width:900px}")),
        shiny::tags$head(shiny::tags$style(".modal-dialog.modal-lg{width:1400px}")),
        shiny::tags$head(shiny::tags$style(".modal-dialog.modal-sm{width:400px}"))
        ##shinyalert::useShinyalert()
    )    
    server = function(input, output, session)
    {
        return_obj <- moduleSERVER(id, ...)
        ## return_obj <- moduleSERVER(id)
        shiny::showModal( shiny::modalDialog(
            moduleUI(id),
            footer = shiny::tagList(
                ## shiny::modalButton("Cancel"),
                shiny::actionButton("gdgt_close","X")
            ),
            size = size,
            easyClose = FALSE,
            fade = FALSE
        ))
        shiny::observeEvent( input$gdgt_close, {
            cat("[gadgetize2] *** closing gadget ***\n")
            shiny::stopApp(return_obj())
        })
    }
    
    pgx <- shiny::runGadget(ui, server)
    ## shiny::shinyApp(ui, server)
    cat(names(pgx))
    pgx
}

this.style <- function(id, css, ns=NULL) {
    if(!is.null(ns)) id <- ns(id)
    shiny::tags$head(shiny::tags$style(paste0("#",id," ",css)))
}

alertDataLoaded <- function(session, ngs) {
    if(!is.null(ngs)) return()
    if(0) {
        shinyWidgets::sendSweetAlert(
                          session = session,
                          ##title = "No dataset loaded",
                          title = NULL,
                          text = "Please first load a dataset"
                      )
    }
    message("[alertDataLoaded] WARNING:: no PGX object")
}

pgx.randomCartoon <- function() {
    cartoon_list <- list(
        list(slogan="Visual analytics. See and understand", img="data-graph-wisdom.jpg"),
        list(slogan="Fasten your seat belts. Accelerated discovery", img="cartoon-speedup.jpg"),
        list(slogan="Analytics anywhere. Anytime.", img="cartoon-cloudservice.jpg"),
        ## list(slogan="Do-it-yourself. Yes you can.", img="bigomics-rockstar3.jpg"),
        list(slogan="Analyze with confidence. Be a rockstar", img="bigomics-rockstar3.jpg"),
        list(slogan="Fast track your Bioinformatics", img="selfservice-checkout2.png"),
        list(slogan="Integrate more. Dig deeper", img="cartoon-integration.jpg"),
        list(slogan="Your analysis doesn't take coffee breaks", img="gone-for-coffee.png"),
        list(slogan="Too much data? Help yourself", img="cartoon-datahelp2.jpg"),    
        list(slogan="Big Friendly Omics", img="big-friendly-omics1.jpg"),
        list(slogan="Big Data meets Biology", img="bigdata-meets.png")
    )
    ##randomCartoon <- shiny::reactive({
    ##invalidateLater(20000)
    cartoon <- sample(cartoon_list,1)[[1]]
    cartoon$img2 = file.path("cartoons",cartoon$img)
    cartoon$img  = file.path("www/cartoons",cartoon$img)
    cartoon
}

pgx.showCartoonModal <- function(msg="Loading data...", img.path="www/cartoons")
{    
    cartoon_list <- list(
        list(slogan="Visual analytics. See and understand", img="data-graph-wisdom.jpg"),
        list(slogan="Fasten your seat belts. Accelerated discovery", img="cartoon-speedup.jpg"),
        list(slogan="Analytics anywhere. Anytime.", img="cartoon-cloudservice.jpg"),
        ## list(slogan="Do-it-yourself. Yes you can.", img="bigomics-rockstar3.jpg"),
        list(slogan="Analyze with confidence. Be a rockstar", img="bigomics-rockstar3.jpg"),
        list(slogan="Fast track your Bioinformatics", img="selfservice-checkout2.png"),
        list(slogan="Integrate more. Dig deeper", img="cartoon-integration.jpg"),
        list(slogan="Your analysis doesn't take coffee breaks", img="gone-for-coffee.png"),
        list(slogan="Too much data? Help yourself", img="cartoon-datahelp2.jpg"),    
        list(slogan="Big Friendly Omics", img="big-friendly-omics1.jpg"),
        list(slogan="Big Data meets Biology", img="bigdata-meets.png")
    )
    
    randomCartoon <- function() {
        cartoon <- sample(cartoon_list,1)[[1]]
        ##cartoon$img = file.path("cartoons",cartoon$img)
        cartoon$img  = file.path(img.path,cartoon$img)
        cartoon
    }
    
    toon <- randomCartoon()
    shiny::showModal(shiny::modalDialog(
        ##title = shiny::HTML("<center><h4>Omics Playground</h4></center>"),
        title = shiny::HTML("<center><h2>",toon$slogan,"</h2><h4>with Omics Playground</h4></center>"),
        shiny::fillRow(
            flex=c(1,NA,1), shiny::br(),
            shiny::img(src = base64enc::dataURI(file=toon$img), width="auto", height="300px"),
            shiny::br()
        ),
        footer = shiny::HTML("<center><p>",msg,"  &nbsp; Please wait</p></center>"),
            size="l", easyClose=FALSE, fade=TRUE))

}

pgx.showSmallModal <- function(msg="Please wait...")
{    
    shiny::showModal(shiny::modalDialog(
        ##title = shiny::HTML("<center><h4>Omics Playground</h4></center>"),
        title = NULL,
        shiny::HTML("<br><center><p>",msg,"</p></center>"),
        footer = NULL,
        size="s", easyClose=FALSE, fade=FALSE))
}
