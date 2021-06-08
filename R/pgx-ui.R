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
    require(miniUI)
    id = sub(".*file","gadget",tempfile())  ## random ID
    ui = miniPage(
        tags$head(tags$style(".modal-dialog{width:900px}")),
        tags$head(tags$style(".modal-dialog.modal-lg{width:1400px}")),
        tags$head(tags$style(".modal-dialog.modal-sm{width:400px}")),        
        useShinyalert(),
        gadgetTitleBar(title),
        miniContentPanel(moduleUI(id))
    )
    server = function(input, output, session) {
        return_obj <- moduleSERVER(id, ...)
        ##return_obj <- moduleSERVER(id)    
        observeEvent( input$done, {
            stopApp(return_obj())
        })
    }
    X <- runGadget(ui, server)
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
    require(shiny)
    require(shinyalert)
    id = sub(".*file","gadget",tempfile())  ## random ID
    ui = fluidPage(
        tags$head(tags$style(".modal-dialog{width:900px}")),
        tags$head(tags$style(".modal-dialog.modal-lg{width:1400px}")),
        tags$head(tags$style(".modal-dialog.modal-sm{width:400px}")),
        useShinyalert()
    )    
    server = function(input, output, session)
    {
        return_obj <- moduleSERVER(id, ...)
        ## return_obj <- moduleSERVER(id)
        showModal( modalDialog(
            moduleUI(id),
            footer = tagList(
                ## modalButton("Cancel"),
                actionButton("gdgt_close","X")
            ),
            size = size,
            easyClose = FALSE,
            fade = FALSE
        ))
        observeEvent( input$gdgt_close, {
            cat("[gadgetize2] *** closing gadget ***\n")
            stopApp(return_obj())
        })
    }
    
    pgx <- runGadget(ui, server)
    ## shinyApp(ui, server)
    cat(names(pgx))
    pgx
}

this.style <- function(id, css, ns=NULL) {
    if(!is.null(ns)) id <- ns(id)
    tags$head(tags$style(paste0("#",id," ",css)))
}


alertDataLoaded <- function(session, ngs) {
    if(!is.null(ngs)) return()
    sendSweetAlert(
        session = session,
        ##title = "No dataset loaded",
        title = NULL,
        text = "Please first load a dataset"
    )
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
    ##randomCartoon <- reactive({
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
    showModal(modalDialog(
        ##title = HTML("<center><h4>Omics Playground</h4></center>"),
        title = HTML("<center><h2>",toon$slogan,"</h2><h4>with Omics Playground</h4></center>"),
        fillRow(
            flex=c(1,NA,1), br(),
            img(src = base64enc::dataURI(file=toon$img), width="auto", height="300px"),
            br()
        ),
        footer = HTML("<center><p>",msg,"  &nbsp; Please wait</p></center>"),
            size="m", easyClose=FALSE, fade=TRUE))

}

pgx.showSmallModal <- function(msg="Please wait...")
{    
    showModal(modalDialog(
        ##title = HTML("<center><h4>Omics Playground</h4></center>"),
        title = NULL,
        HTML("<br><center><p>",msg,"</p></center>"),
        footer = NULL,
        size="s", easyClose=FALSE, fade=FALSE))
}
