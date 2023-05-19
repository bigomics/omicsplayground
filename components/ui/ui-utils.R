##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#star.symbols <- function(n, pch="\u2605") {
#    if(n==0) return("")
#    paste(rep(pch,n),collapse="")
#}

addWatermark.PDF <- function(file) {
    if(system("which pdftk",ignore.stdout=TRUE)==1) return ## if no pdftk installed...
    mark <- file.path(FILES,"watermark.pdf")
    tmp <- paste0(gsub("file","plot",tempfile()),".pdf")
    cmd <- paste("pdftk",file,"stamp",mark,"output",tmp) ## NEED pdftk installed!!!
    cmd
    system(cmd)
    file.copy(tmp,file,overwrite=TRUE)
    unlink(tmp)
}

addWatermark.PNG2 <- function(file,
                              w, h, mt=0.09, out=file,
                              mark=file.path(FILES,"watermark-logo.png"),
                              logo.scale=0.045) {
    if(system("which convert",ignore.stdout=TRUE)==1) return ## if no pdftk installed...
    tmp <- paste0(gsub("file","plot",tempfile()),".png")
    logo.height = max(logo.scale*h , logo.scale*w/3)
    logo.width = logo.height * 4
    logo.x = 0.2 * logo.width
    logo.y = 0.3 * logo.height
    top.margin = mt*h
    if(top.margin>0) {
      cmd.str = "convert %s -extent %.0fx%.0f+%.0f-%.0f %s"
      cmd0 = sprintf( cmd.str, file, w, h+top.margin, 0, top.margin, file)
      system(cmd0)
    }
    cmd.str = "convert %s \\( %s -thumbnail %.0fx%.0f \\) -geometry +%.0f+%.0f -composite %s"
    cmd = sprintf( cmd.str, file, mark, logo.width, logo.height, logo.x, logo.y, tmp)
    system(cmd)
    file.copy(tmp,out,overwrite=TRUE)
    unlink(tmp)
}

addWatermark.PDF2 <- function(file, w, h, out=file, mt=0.09,
                              mark=file.path(FILES,"watermark-logo.pdf"),
                              logo.scale=0.045 ) {
    if(system("which pdftk",ignore.stdout=TRUE)==1) return ## if no pdftk installed...    
    tmp1 <- paste0(gsub("file","plot",tempfile()),".pdf")
    tmp2 <- paste0(gsub("file","plot",tempfile()),".pdf")    
    tmp3 <- paste0(gsub("file","plot",tempfile()),".pdf")    
    logo.height = max(logo.scale*h , logo.scale*w/3) * 720    
    logo.width = logo.height * 4
    scale.cmd = sprintf("gs -o %s -sDEVICE=pdfwrite  -dDEVICEWIDTH=%0.f -dDEVICEHEIGHT=%0.f -dPDFFitPage -dAutoRotatePages=/None -f %s", tmp1, logo.width, logo.height, mark)
    system(scale.cmd, ignore.stdout=TRUE, ignore.stderr=FALSE)
    ppi = 720
    if(mt>0) {
      ## add margin??      
      h = h + h*mt
      cmd0 = sprintf("gs -o %s -sDEVICE=pdfwrite -g%.0fx%.0f -dPDFFitPage -f -dAutoRotatePages=/None %s",tmp3, w*ppi, 1.03*h*ppi, file)      
      system(cmd0, ignore.stdout=FALSE, ignore.stderr=FALSE)
    } else {
      file.copy(file,tmp3,overwrite=TRUE)
    }
    ## create empty page with same plot size and logo translated to topleft
    cmd1 = sprintf("gs -o %s -sDEVICE=pdfwrite -g%.0fx%.0f -c '<</PageOffset [%.0f %.0f]>> setpagedevice' -f %s",
      tmp2, w*720, h*720, 0.15*logo.width/10, h*72-1*logo.height/10, tmp1)
    system(cmd1, ignore.stdout=TRUE, ignore.stderr=FALSE)

    ## merge: overlay logo and plot
    cmd3 <- paste("pdftk",tmp3,"stamp",tmp2,"output",out) ## NEED pdftk installed!!!    
    system(cmd3)
   ## file.copy(tmp3,file,overwrite=TRUE)
    unlink(tmp1)
    unlink(tmp2)
    unlink(tmp3)    
}

gadgetize <- function(moduleUI, moduleSERVER, title="shiny gadget", ...)
{
    ## Creates gadget from a Shiny module. Gadget are browser-based UI
    ## applets for a single task that can be run from the R command
    ## line. It is not to be used inside Shiny programs themselves.
    ##

    id = sub(".*file","gadget",tempfile())  ## random ID
    ui = miniUI::miniPage(
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
        cartoon$img2 = paste0("static/cartoons/",cartoon$img)
        cartoon$img  = file.path(img.path,cartoon$img)
        cartoon
    }

    toon <- randomCartoon()
    shiny::showModal(shiny::modalDialog(
        #title = shiny::HTML("<center><h4>Omics Playground</h4></center>"),
        ##title = shiny::HTML("<center><h2>",toon$slogan,"</h2><h4>with Omics Playground</h4></center>"),
        title = shiny::div( shiny::h2(toon$slogan), shiny::p("with Omics Playground"), style="text-align:center;"),
        #title = div(h2(toon$slogan), style="text-align:center;width:100%;"),
        shiny::img(src = toon$img2, class = "img-fluid"),
        ##footer = div(msg, style="text-align:center;"),
        footer = fillRow( flex=c(1,NA,1), " ", msg, " "),
        size = "l",
        easyClose = FALSE,
        fade = TRUE
    ))
    ## Sys.sleep(600)  # for debugging
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
