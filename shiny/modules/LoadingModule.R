LoadingInputs <- function(id) {
    ns <- NS(id)  ## namespace
    tagList(
        uiOutput(ns("description"))
        ## uiOutput(ns("inputsUI"))
        ## uiOutput(ns("socialButtons"))
    )
}

LoadingUI <- function(id) {
    ns <- NS(id)  ## namespace
    fillCol(
        height = 750,
        flex = c(NA,1),
        fillCol(
            height=130,
            uiOutput(ns("valueboxes_UI"))
        ),
        tabsetPanel(
            id = ns("tabs"),
            tabPanel("Public datasets",uiOutput(ns("pgxtable_UI"))),
            tabPanel("Upload data",uiOutput(ns("upload_UI"))),
            tabPanel("Visitors map",uiOutput(ns("usersmap_UI"))),
            tabPanel("Community forum",uiOutput(ns("forum_UI")))
        )
    )
}


LoadingModule <- function(input, output, session, hideModeButton=TRUE,
                          max.limits=c("samples"=1000,"comparisons"=20,"genes"=19999),
                          defaultMode="BASIC", authentication="none")
{
    ns <- session$ns ## NAMESPACE

    message("[LoadingModule] DEBUG=",DEBUG)
    message("[LoadingModule] USER_MODE=",USER_MODE)
    
    ##useShinyjs(rmd=TRUE)
    useShinyjs()
    useSweetAlert()
    SHOWSPLASH=TRUE
    ## SHOWSPLASH=FALSE
    hideModeButton <- toupper(hideModeButton)

    LOGIN_AUTHENTICATION = authentication
    ## LOGIN_AUTHENTICATION = "none"
    ##LOGIN_AUTHENTICATION = "register"
    ##LOGIN_AUTHENTICATION = "password"
    CREDENTIALS = NULL
    if(file.exists("CREDENTIALS")) {
        CREDENTIALS <- read.table("CREDENTIALS",row.names=1,
                                  header=TRUE, stringsAsFactors=FALSE)
        LOGIN_AUTHENTICATION = "password"
    }

    auth <- callModule(
        AuthenticationDialog, "auth",
        type = LOGIN_AUTHENTICATION,
        credentials = CREDENTIALS
    )

    ##-----------------------------------------------------------------------------
    ## Description
    ##-----------------------------------------------------------------------------
    description = "<b>Omics Playground</b> is an interactive self-service bioinformatics platform for the analysis, visualization and interpretation of transcriptomics and proteomics data. Life scientists can easily perform complex data analysis and visualization without coding, and significantly reduce the time to discovery."
    
    output$description <- renderUI(HTML(description))

    ##-----------------------------------------------------------------------------
    ## Show current dataset on each page
    ##-----------------------------------------------------------------------------
    curDataSet <- reactive({
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)
        ##HTML("<b>dataset :</b>",ngs$name,"")
        ##HTML("<h3>",ngs$name,"</h3>")
        dname <- gsub("^.*/|[.]pgx","",ngs$name)
        HTML("<div class='current-data'>",dname,"</div>")
    })
    output$current_dataset <- renderText({ curDataSet() })

    ##-----------------------------------------------------------------------------
    ## User interface
    ##-----------------------------------------------------------------------------

    usermode_infotip = "Select BASIC or PRO user mode. The BASIC mode should be sufficient for most users. The PRO mode unlocks more algorithms, extra visualization panels and advanced analysis modules for single-cell RNA-seq, cell profiling and biomarker analysis."

    USERLEVELS = c("BASIC","PRO")
    ## if(DEV.VERSION) USERLEVELS = c("BASIC","PRO","DEV")
    ##USERMODE <- reactiveVal( factor("BASIC",levels=USERLEVELS) )
    USERMODE <- reactiveVal( factor(toupper(defaultMode),levels=USERLEVELS) )

    output$inputsUI <- renderUI({

        usermodeUI <- tipify(radioGroupButtons(
            inputId = ns("main_usermode"),
            label = "User mode:",
            choices = USERLEVELS,
            selected = toupper(defaultMode),
            status = "warning",
            checkIcon = list(yes = icon("ok", lib = "glyphicon"))),
            usermode_infotip, placement="bottom", options = list(container="body")
            )
        
        if(hideModeButton) usermodeUI <- NULL
        
        ui <- tagList(
            usermodeUI,
            ##br(),br(),br(),
            p(strong("Dataset info:")),
            div( htmlOutput(ns("dataset_info")), id="datainfo"),
            br(),
            tipify( actionButton(ns("loadbutton"),label="Load dataset",class="load-button"),
                   "Click to load the selected dataset.", placement="bottom")
        )
        ui
    })
    outputOptions(output, "inputsUI", suspendWhenHidden=FALSE)
    
    output$dataset_info <- renderText({
        sec <- currentSection()
        inf <- selectedDataSetInfo()
        inf["conditions"] <- gsub("[,]"," ",inf["conditions"])
        inf <- sapply(inf, function(s) substring(s,1,500)) 
        if(sec=="upload-data") {
            HTML(paste("<p>Please upload dataset<br>"))
        } else if(length(inf)==0) {
            HTML(paste("<p>Please select a dataset<br>"))
        } else {
            HTML(paste("<p><b>",names(inf),"</b>:", inf,"<br>"))
        }
    })

    ##-----------------------------------------------------------------------------
    ## UPDATING PGX FILE INFO
    ##-----------------------------------------------------------------------------
    ##PGXINFO <- pgx.updateInfoFile(PGX.DIR, file="datasets-info.csv", 
    ##                           force=FALSE, verbose=TRUE )
    PGXINFO <- pgx.scanInfoFile(PGX.DIR, file="datasets-info.csv", verbose=TRUE )
    dim(PGXINFO)
    
    ##=================================================================================
    ##========================== MODAL DIALOGS ========================================
    ##=================================================================================

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
        ##randomCartoon <- reactive({
        ##invalidateLater(20000)
        cartoon <- sample(cartoon_list,1)[[1]]
        cartoon$img2 = file.path("cartoons",cartoon$img)
        cartoon$img  = file.path("www/cartoons",cartoon$img)
        cartoon
    }##)
    
    startup_count=0
    ##require(shinyparticles)
    require(particlesjs)
    particlesjs.conf <- rjson::fromJSON(file="resources/particlesjs-config.json")

    showStartupModal <- function(once=FALSE) {    

        if(length(input$loadbutton)==0) {
            dbg("[showStartupModal] UI not ready. skipping")
            ##delay(8000, selectRows(proxy = dataTableProxy("pgxtable"), selected=1))
            ##delay(8000, shinyjs::click("loadbutton"))
            return(NULL)  ## UI not ready???
        }
        if(once && startup_count>0) return(NULL)
        
        dbg("showStartupModal: showing!\n")    
        showModal(modalDialog(
            id="modal-splash",
            div(
                id = "particles-target",
                img(src = base64enc::dataURI(file="www/splash.png"),
                    width="100%", height="auto%", style = "position:absolute;"),
                ##div("Big Omics Data", class="splash-title"),
                ##div("Isn't big anymore with Omics Playground", class="splash-subtitle"),
                ##style="height: 500px; width: 100%; background-color: #2c81e2;"
                style="height: 500px; width: 100%;"                
                ##height="500px", style="margin-left: -14px; margin-top: -4px;"
            ),            
            footer = tagList(
                actionButton("action1","Read-the-docs", icon=icon("book"),
                             onclick="window.open('https://omicsplayground.readthedocs.io','_blank')"),
                actionButton("action2","Watch the tutorial", icon=icon("youtube"),
                             onclick="window.open('https://www.youtube.com/watch?v=_Q2LJmb2ihU&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-','_blank')"),
                actionButton("action3","Get the source", icon=icon("github"),
                             onclick="window.open('https://github.com/bigomics/omicsplayground','_blank')"),
                actionButton("action4","Docker image", icon=icon("docker"),
                             onclick="window.open('https://hub.docker.com/r/bigomics/omicsplayground','_blank')"),
                actionButton("action5","User forum", icon=icon("users"),
                             onclick="window.open('https://groups.google.com/d/forum/omicsplayground','_blank')"),
                ## actionButton(ns("action_beer"),"Buy us a beer!", icon=icon("beer")),
                ## actionButton(ns("action_pizza"),"Donate pizza", icon=icon("pizza-slice")),
                ## modalButton("Let's start!")
                actionButton(ns("action_play"), "Let's play!")
            ),
            particles(config=particlesjs.conf, target_id = "particles-target", timeout = 1000),
            size="l", easyClose=FALSE, fade=FALSE))
        ##startup_count <<- startup_count + 1    
        dbg("showStartupModal done!\n")
    }

    ## observeEvent( input$action_beer, {
    ##     dbg("buy beer button action\n")
    ##     startup_count <<- startup_count + 1    
    ##     USER$logged <- TRUE
    ##     USER$name   <- "beer buddy"
    ##     removeModal()
    ##     ##alert("Wow. Thanks buddy!")
    ##     sendSweetAlert(
    ##         session=session, title="Wow. Thanks buddy!",
    ##         text = "Free entrance for you!", type = "info")
    ##     ##Sys.sleep(4);removeModal()
    ## })

    observeEvent( input$action_play, {

        dbg("action_play:: play button action\n")
        cat("action_play:: LOGIN_AUTHENTICATION=",LOGIN_AUTHENTICATION,"\n")

        startup_count <<- startup_count + 1    
        if(LOGIN_AUTHENTICATION!="none") {
            message("[LoadingModule] authentication required")
            ## showLogin()  ## $
            AuthenticationUI(ns("auth"))
        } else {
            message("[LoadingModule] continuing without authentication")
            removeModal()
        }    

    })

    selectedDataSetInfo <- reactive({
        ##sel <- input$pgxtable_rows_selected
        sel <- pgxtable$rows_selected()
        if(is.null(sel) || length(sel)==0) return(NULL)
        df <- getPGXTable()
        unlist(lapply(df[sel,],as.character))
    })

    selectedDataSet <- reactive({
        ##sel <- input$pgxtable_rows_selected
        sel <- pgxtable$rows_selected()
        if(is.null(sel) || length(sel)==0) return(NULL)
        df <- getPGXTable()
        as.character(df$dataset[sel])
    })

    output$loading_image <- renderImage({
        toon <- randomCartoon()
        list(src = toon$img,
             contentType = 'image/png',
             width = "100%", height = "100%", ## actual size: 1040x800         
             alt = "loading image")
    }, deleteFile=FALSE)        


    showloading_ntime=0
    showLoadingModal <- function(msg="Loading data...") {
        toon <- randomCartoon()
        showModal(modalDialog(
            ##title = HTML("<center><h4>Omics Playground</h4></center>"),
            title = HTML("<center><h2>",toon$slogan,"</h2><h4>with Omics Playground</h4></center>"),
            fillRow(flex=c(1,NA,1), br(),
                    ##img0,
                    ##imageOutput("loading_image", width="auto", height="250px"),
                    img(src = base64enc::dataURI(file=toon$img), width="auto", height="300px"),
                    br()),
            footer = HTML("<center><p>",msg,"  &nbsp; Please wait</p></center>"),
            size="l", easyClose=FALSE, fade=TRUE))
        ## Sys.sleep(5)
        showloading_ntime <<- 1
        dbg("showLoadingModal done!\n")
    }

    currentSection <- reactive({
        cdata <- session$clientData
        sub("section-","",cdata[["url_hash"]])
    })


    ##=================================================================================
    ##==================== USER AUTHENTICATION ========================================
    ##=================================================================================

    ## USER <- reactiveValues( logged = FALSE, name="anonymous")
    
    ##=================================================================================
    ##======================== USER LEVEL =============================================
    ##=================================================================================

    ##-----------------------------------------------------------------------------
    ## Select UI user level
    ##-----------------------------------------------------------------------------

    observeEvent( input$main_usermode,
    {
        ## Observe the usermode button and switch levels
        ##
        ##
        usermode <- input$main_usermode
        if(is.null(usermode)) usermode <- "BASIC"

        cat(">>> switching USERMODE to",usermode,"\n")
        
        if(DEV.VERSION) {            
            ##dbg("VALID OUTPUT.OPTIONS =",names(outputOptions(output)))
            ##dbg("VALID OUTPUT NAMES =",names(output))
            dbg("VALID INPUT NAMES =",names(input))
            dbg("session$clientData NAMES =",names(session$clientData))
        }
        
        if(usermode == "BASIC") {
            dbg("observeEvent::main_usermode : switching to BASIC mode")
            ##shinyjs::html("navbar-brand","Omics Playground (basic)")
            shinyjs::hide(selector = "div.download-button")
            shinyjs::hide(selector = "div.modebar")
            shinyjs::hide(selector = "div.pro-feature")
            USERMODE("BASIC")
        }
        
        if(usermode %in% c("PRO","DEV")) {
            dbg("observeEvent::main_usermode : switching to PRO mode")
            shinyjs::show(selector = "div.download-button")
            shinyjs::show(selector = "div.modebar")
            shinyjs::show(selector = "div.pro-feature")
            USERMODE("PRO")
        }

        if(usermode %in% c("DEV")) {
            dbg("observeEvent::main_usermode : switching to DEV mode")
            ##shinyjs::show(selector = "div.download-button")
            USERMODE("DEV")
        }
        
        ##updateNavbarPage(session, "navbar", selected=1)
    }, ignoreNULL=TRUE, ignoreInit=TRUE )

    ##================================================================================
    ##====================== INPUT DATA REACTIVE OBJECT ==============================
    ##================================================================================

    currentPGX <- reactiveVal(NULL)

    ##inputData <- eventReactive( reload(), {
    inputData <- reactive({
        ##-----------------------------------------------------------------
        ## This is the main loader function that loads the ngs object.
        ##-----------------------------------------------------------------
        dbg("[LoadingModule::inputData] ---------- reacted ---------------\n")
        dbg("[LoadingModule::inputData] LOGIN_AUTHENTICATION=",LOGIN_AUTHENTICATION,"\n")
        dbg("[LoadingModule::inputData] auth$logged=",auth$logged(),"\n")
        
        ## authenicate user if needed
        ## if(!USER$logged) showLogin()
        ## if(LOGIN_AUTHENTICATION!="none" && USER$logged) showLogin()
        if(LOGIN_AUTHENTICATION!="none" && !auth$logged()) return(NULL)

        pgx <- currentPGX()
        dbg("[LoadingModule::inputData] is.null(pgx)=",is.null(pgx),"\n")        
        return(pgx)
    })

    observeEvent( input$loadbutton, {

        ## Observe button press
        btn <- isolate(input$loadbutton)
        pgx = NULL
        pgx = isolate(selectedDataSet())

        dbg("[LoadingModule::<loadbutton>] loadbutton=",btn,"\n")
        dbg("[LoadingModule::<loadbutton>] pgx.selected=",pgx)

        if(!is.null(btn) && btn!=0 && !is.null(pgx)) {
            ## show loading pop-up
            showLoadingModal()
        }

        if(is.null(pgx) || pgx=="" || length(pgx)==0) {
            ## Set to default data set if table is not ready yet.
            ## pgx <- PGXINFO$dataset[1]
            ## dbg("[observe:loadbutton] setting to default data set = ",pgx,"\n")
            return(NULL)
        }

        
        pgx.path <- PGX.DIR[file.exists(file.path(PGX.DIR,pgx))]
        pgx1 = file.path(pgx.path,pgx)
        pgx1
        if(file.exists(pgx1)) {
            dbg("[LoadingModule::<loadbutton>] LOADING",pgx1,"\n")
            ##withProgress(message='loading...', value=0.8,
            load(pgx1,verbose=0)
        } else {
            cat("[LoadingModule::<loadbutton>] ERROR file not found : ",pgx1,"\n")
            removeModal()
            return(NULL)
        }
        
        ##----------------- update input
        dbg("[LoadingModule::<loadbutton>] head.names.NGS=",head(names(ngs)))
        dbg("[LoadingModule::<loadbutton>] initializing PGX object")
        ngs <- pgx.initialize(ngs)
        if(is.null(ngs)) {
            cat("[LoadingModule::<loadbutton>] ERROR in object initialization\n")
            showNotification("ERROR in object initialization!\n")
            removeModal()
            return(NULL)
        }
        if(is.null(ngs$name)) ngs$name <- sub("[.]pgx$","",pgx)

        ##----------------- remove modal??
        if(startup_count>0) {
            Sys.sleep(4)
            removeModal()
        }
        
        currentPGX(ngs)
        dbg("[LoadingModule::<loadbutton>] ready! \n")
    })
    ##}, ignoreNULL=FALSE )
    ##}, ignoreNULL=TRUE )


    ## ================================================================================
    ## ===================== VALUE BOXES UI ===========================================
    ## ================================================================================

    require(shinydashboard)
    ## useShinydashboard()
    vbox <- function(value, label) {
        box(
            h1(value, style="font-weight: 800; color: white; padding: 16px 0 0 0; margin: 0 0 0 20px;"),
            h5(label, style="margin: 0 0 0 20px; font-weight: 400; color: white; padding-bottom: 25px"),
            ##h1(value, style="font-weight: 800; padding: 16px 0 0 0; margin: 0 0 0 20px;"),
            ##h5(label, style="margin: 0 0 0 20px; font-weight: 400; padding-bottom: 25px"),
            width="100%", class="vbox")
    }

    output$valuebox1 <- renderUI({
        pgx <- getPGXTable()
        req(pgx)
        ndatasets = "..."
        ndatasets <- nrow(pgx)
        vbox( ndatasets, "data sets")     
    })
    output$valuebox2 <- renderUI({
        pgx <- getPGXTable()
        req(pgx)
        ##dbg("valuebox2:: pgx$nsamples=",pgx$nsamples)
        nsamples <- sum(as.integer(pgx$nsamples),na.rm=TRUE)
        vbox( nsamples, "number of samples")     
    })
    output$valuebox3 <- renderUI({
        pgx <- getPGXTable()
        req(pgx)
        ##dbg("valuebox3:: pgx$nsamples=",pgx$nsamples)
        nvalues <- sum(as.integer(pgx$nsamples) * (as.integer(pgx$ngenes)
            + as.integer(pgx$nsets)),na.rm=TRUE)
        nvalues1 <- format(nvalues, nsmall=, big.mark=" ")  # 1,000.6
        vbox( nvalues1, "data points") 
    })

    output$valueboxes_UI <- renderUI({
        fillRow(
            height=115,
            uiOutput(ns("valuebox1")),
            uiOutput(ns("valuebox2")), 
            uiOutput(ns("valuebox3"))
        )
    })

    ##================================================================================
    ## Data sets
    ##================================================================================
    
    ##split=" ";n=5
    andothers <- function(s, split=" ", n=8) {
        if(is.na(s)) return("")
        s <- sub("^[ ]*","",s)
        s <- sub("[ ]+"," ",s)
        s1 <- strsplit(s, split=split)[[1]]
        if(length(s1)<=n) return(s)
        n2 <- setdiff(length(s1),n)
        paste(paste(head(s1,n), collapse=" "),"(+",n2,"others)")
    }

    getPGXTable <- reactive({
        ## get table of data sets
        ##
        ##
        
        if(is.null(PGXINFO)) return(NULL)    
        dbg("[LoadingModule:getPGXTable] reacted")
        
        df <- PGXINFO
        pgx.files = dir(PGX.DIR, pattern=".pgx$")
        sel <- sub("[.]pgx$","",df$dataset) %in% sub("[.]pgx$","",pgx.files)
        df <- df[sel,,drop=FALSE]
        
        ##kk = unique(c("dataset","datatype","organism","description",colnames(df)))
        kk = unique(c("dataset","datatype","organism","description","nsamples",
                      "ngenes","nsets","conditions","date"))
        kk = intersect(kk,colnames(df))
        df = df[,kk]               
        df = df[order(df$dataset),]   ## sort alphabetically...
        rownames(df) <- NULL
        df
    })


    pgxTable.RENDER <- reactive({

        if(SHOWSPLASH) showStartupModal(once=TRUE)
        
        dbg("<pgxTable.RENDER> reacted")

        df <- getPGXTable()
        req(df)

        df$dataset  <- gsub("[.]pgx$"," ",df$dataset)
        df$conditions  <- gsub("[,]"," ",df$conditions)
        df$conditions  <- sapply(as.character(df$conditions), andothers, split=" ", n=5)
        df$description <- shortstring(as.character(df$description),200)
        df$nsets <- NULL
        df$date  <- NULL
        
        DT::datatable( df,
                      class = 'compact cell-border stripe hover',
                      rownames=TRUE,
                      extensions = c('Scroller'),
                      selection = list(mode='single', target='row', selected=NULL ),
                      ## filter = "top",
                      fillContainer = TRUE,
                      options=list(
                          ##dom = 'Blfrtip',
                          dom = 'frti',
                          ##columnDefs = list(list(searchable = FALSE, targets = 1)),
                          pageLength = 1000, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = FALSE,
                          ##scrollY =400, ## scroller=TRUE,
                          scrollY = '100vh', ## scroller=TRUE,
                          deferRender=TRUE
                      )  ## end of options.list 
                      )  %>%
            DT::formatStyle(0, target='row', fontSize='11.5px', lineHeight='95%')

    })

    pgxDatasetOverview.RENDER <- reactive({        
        df <- getPGXTable()
        req(df)
        ##sel = input$pgxtable_rows_selected
        sel <- pgxtable$rows_selected()
        sel.dataset = df$dataset[sel]
        par(mfrow=c(1,3))
        plot(sin)
        plot(cos)
        plot(exp)        
    })

    pgxtable_text = "This table contains a general information about all available datasets within the platform. For each dataset, it reports a brief description as well as the total number of samples, genes, gene sets (or pathways), corresponding phenotypes and the collection date."

    ## pgxtable_module <- tableModule(
    ##     id = "pgxtable",
    ##     func = pgxTable.RENDER,
    ##     info.text = pgxtable_text,
    ##     title="Datasets"
    ## )    
    ## output <- attachModule(output, pgxtable_module)
    pgxtable <- callModule(
        tableModule, id = "pgxtable",
        func = pgxTable.RENDER,
        title = "Datasets",
        height = 580)
    ##outputOptions(output, "pgxtable", suspendWhenHidden=FALSE) ## important!

    output$pgxtable_UI <- renderUI({    
        fillCol(
            height = 750,
            fillRow(
                flex = c(1,0.1,4.5),
                wellPanel(
                    uiOutput(ns("inputsUI"))
                ),
                br(),
                tableWidget(ns("pgxtable"))
            )
        )        
    })

    ##================================================================================
    ## Upload data
    ##================================================================================

    ## Some 'global' reactive variables used in this file
    uploaded_files <- reactiveValues()
    ##shinyjs::disable("upload_compute")
    
    output$downloadExampleData <- downloadHandler(
        filename = "exampledata.zip",
        content = function(file) {
            zip = file.path(FILES,"exampledata.zip")
            file.copy(zip,file)
        }
    )

    upload_info = "<h4>User file upload</h4><p>Please prepare the data files in CSV format as listed below. It is important to name the files exactly as shown. The file format must be comma-separated-values (CSV) text. Be sure the dimensions, rownames and column names match for all files. You can download a zip file with example files here: EXAMPLEZIP. You can upload a maximum of <u>LIMITS</u>. After uploading, in sidebar on the left, provide a name for your dataset. Finally, hit the compute button. The computations may take 10 to 30 minutes depending on the size of your dataset."
    DLlink = downloadLink(ns("downloadExampleData"),"exampledata.zip")
    upload_info = sub("EXAMPLEZIP", DLlink, upload_info)
    
    upload_info2 =
        "<br><h4>Uploaded datasets</h4><p>Below are your uploaded datasets. As a free user, you can only have a maximum of one private dataset. If you want to analyze a new dataset, you must either delete your old dataset or make the dataset public.<br><br>"

    ##upload_filetypes = c("text/csv","text/comma-separated-values,text/plain",".csv")
    upload_filetypes = c(".csv",".pgx")    
    output$upload_UI <- renderUI({    

        limits <- paste(max.limits["samples"],"samples and",
                        max.limits["comparisons"],"comparisons")
        
        upload_info = sub("LIMITS", limits, upload_info)
        ##if(USERMODE()=="PRO") upload_info = sub("LIMITS", pro.limits, upload_info)    
        userdataUI <- NULL
        if(DEV.VERSION) {
            userdataUI <- DT::dataTableOutput(ns("userDatasetsUI"))
        }
        
        fillCol(
            flex = c(1,1.5),
            height = 750,
            fillRow(
                flex = c(1,0.1,4.5),
                wellPanel(
                    fileInput(ns("upload_files"), "Choose files",
                              multiple = TRUE, accept = upload_filetypes),
                    ## textInput(ns("upload_name"),"Name of dataset:"),
                    ## textAreaInput("upload_description", "Description:", value = NULL,
                    ##              rows=5, placeholder="Describe your data set (minimum 100 characters)"),
                    actionButton(ns("upload_compute"),"Compute!",icon=icon("running"),
                                 class="run-button")
                ),br(),
                fillCol(
                    flex = c(NA,1),
                    div(HTML(upload_info),style="font-size: 13px;"),
                    tableOutput(ns("uploadStatusTableOutput"))
                )
            ),
            fillRow(
                userdataUI
            )
        )
    })

    observeEvent( input$upload_compute, {

        ## are you sure? Any message/warning before computation is done.
        has.pgx <- ("uploaded.pgx" %in% names(uploaded_files))
        has.csv <- all(c("counts.csv","samples.csv","contrasts.csv") %in% names(uploaded_files))
        if(has.pgx) has.pgx <- has.pgx && !is.null(uploaded_files[["uploaded.pgx"]])

        dbg("[observeEvent::upload_compute] names(uploaded_files)=",names(uploaded_files))
        dbg("[observeEvent::upload_compute] has.pgx=",has.pgx)
        dbg("[observeEvent::upload_compute] has.csv=",has.csv)
        
        if(!has.pgx && !has.csv ) {
            message("[LoadingModule::*upload_compute} WARNING: ***** no PGX, no CSV files *****")
            ##removeModal()
            return(NULL)
        }

        showModal( modalDialog(
            HTML("Your data will be uploaded for computation. By uploading the data you accept our EULA and you confirm that the data has been anonymized. For research use only."),
            title = NULL,
            size = "s",
            footer = tagList(
                modalButton("Cancel"),
                actionButton(ns("myconfirmation"),"Confirm", icon=NULL)
            )
        ))
        
    })

    observeEvent( input$myconfirmation, {
        ## 
        ## Start pre-computing the object from the uploaded files
        ## after confirmation is received.
        ##
        
        dbg("[LoadingModule::*myconfirmation] names(uploaded_files)=",names(uploaded_files))        
        
        ## --------------------- OK start ---------------------------
        has.pgx <- ("uploaded.pgx" %in% names(uploaded_files))       
        if(has.pgx) has.pgx <- has.pgx && !is.null(uploaded_files[["uploaded.pgx"]])

        if(has.pgx) {            
            message("[LoadingModule] ***** using uploaded PGX ******")        
            ngs <- uploaded_files[["uploaded.pgx"]]
            currentPGX(ngs)  ## copy to global reactive variable
        } else {

            message("[LoadingModule] ***** computing from CSV files *****")
            showLoadingModal("Calculating... it's a good time to get a coffee now.")
            
            counts    <- as.matrix(uploaded_files[["counts.csv"]])
            samples   <- data.frame(uploaded_files[["samples.csv"]],stringsAsFactors=FALSE)
            contrasts <- as.matrix(uploaded_files[["contrasts.csv"]])
            contrasts[is.na(contrasts)] <- 0
            
            max.genes = as.integer(max.limits["genes"])
            max.genesets = 9999
            
            if( FALSE && USERMODE() == "BASIC") {
                dbg("[LoadingModule::*myconfirmation] setting BASIC methods")            
                gx.methods   = c("ttest.welch","ttest.rank","trend.limma") ## fastest 3
                gset.methods = c("fisher","gsva","camera")  ## fastest 3            
                ## gx.methods   = c("trend.limma","edger.qlf","edger.lrt")
                ## gset.methods = c("fisher","gsva","fgsea")
                extra.methods = c("meta.go","infer","drugs","wordcloud")
            } else {
                dbg("[LoadingModule::*myconfirmation] setting PRO methods")  
                gx.methods   = c("ttest.welch","trend.limma","edger.qlf","deseq2.wald")
                gset.methods = c("fisher","gsva","fgsea","camera","fry")
                extra.methods = c("meta.go","infer","deconv","drugs-combo",
                                  "wordcloud","connectivity")
                if(ncol(counts) > 750) {
                    ## probably scRNA-seq... to long
                    gx.methods   = c("ttest","ttest.welch","trend.limma") ## only t-test...
                    gset.methods = c("fisher","gsva","fgsea")
                    extra.methods = c("meta.go","infer","deconv","drugs-combo",
                                      "wordcloud","connectivity")
                    max.genes = 10000
                }
            }
            
            ##----------------------------------------------------------------------
            ## Upload and do precomputation
            ##----------------------------------------------------------------------
            start_time <- Sys.time()
            ## Create a Progress object
            progress <- shiny::Progress$new()
            on.exit(progress$close())    
            progress$set(message = "Processing", value = 0)

            progress$inc(0.01, detail = "parsing data")            
            ngs <- pgx.createPGX(
                counts, samples, contrasts, ## genes, 
                only.hugo = TRUE, only.proteincoding = TRUE)
            names(ngs)
            
            ngs <- pgx.computePGX(
                ngs,
                max.genes = max.genes,
                max.genesets = max.genesets, 
                gx.methods = gx.methods,
                gset.methods = gset.methods,
                extra.methods = extra.methods,
                lib.dir = FILES, do.cluster=TRUE,
                progress=progress)
            
            ## ngs <- pgx.computeObjectPGX(
            ##     counts, samples, contrasts,
            ##     max.genes = max.genes,
            ##     gx.methods = gx.methods,
            ##     gset.methods = gset.methods,
            ##     extra.methods = extra.methods,
            ##     lib.dir = FILES, only.hugo = TRUE,
            ##     progress = progress
            ## )

            end_time <- Sys.time()
            delta_time  = end_time - start_time
            delta_time
            dbg("upload_compute :: total processing time of",delta_time,"secs")
            
            names(ngs)
            head(ngs$samples)
            ngs$datatype = "generic"
            ngs$description = "not available"
            ngs.name = "(uploaded)"
            ##ngs.name = gsub("[ ]","-",input$upload_name)
            ngs$name = ngs.name
            ## ngs$datatype = input$upload_datatype
            ## ngs$description = input$upload_description
            ngs$date = date()
        }

        ## initialize and update global PGX object
        ngs <- pgx.initialize(ngs)
        currentPGX(ngs)  ## copy to global reactive variable
        selectRows(proxy = dataTableProxy(ns("pgxtable")), selected=NULL)
        ## shinyjs::click("loadbutton")    

        removeModal()
        showModal( modalDialog(
            HTML("<b>Ready!</b><br>You can now start exploring your data. Tip: to avoid computing again, download your data object locally or save it to the cloud."),
            title = NULL,
            size = "m",
            footer = tagList(
                downloadButton(ns("downloadPGX"), "Download locally", icon=icon("download")),
                actionButton(ns("savedata"), "Save to cloud", icon=icon("save")),
                actionButton(ns("sharedata"), "Share with others", icon=icon("share-alt")),
                modalButton("Start!")
            )
        ))

        ## clean up uploaded_file object
        for(s in names(uploaded_files)) uploaded_files[[s]] <- NULL
                
    })

    output$downloadPGX <- downloadHandler(
        filename = "userdata.pgx",
        content = function(file) {
            ngs <- currentPGX()  ## current dataset
            temp <- tempfile()
            save(ngs, file=temp)
            file.copy(temp,file)
        }
    )

    observeEvent( input$sharedata, {
        dbg("[LoadingModule] observeEvent:input$sharedata reacted")
        showModal( modalDialog(
            HTML("<center><b>New feature!</b><br>Sharing your data with others will be soon available as new feature!</center>"),
            title = NULL, size = "s", fade=FALSE
        ))
        return(NULL)
    })
    
    observeEvent( input$savedata, {

        dbg("[LoadingModule] observeEvent:savedata reacted")
        
        showModal( modalDialog(
            HTML("<center><b>New feature!</b><br>Saving your data to the cloud will be soon available as new feature for Pro users!</center>"),
            title = NULL, size = "s", fade=FALSE
        ))
        return(NULL)
        
        ## -------------------- save PGX file/object
        saving.ok = ( USERMODE() == "PRO")
        if(saving.ok) {

            if(0) {
                ##!!!!!!!!!!!!!!!!!!!!!!!!
                ##!!!!!!! BROKEN !!!!!!!!!
                ##!!!!!!!!!!!!!!!!!!!!!!!!
                ngs.name1 = sub("[.]pgx$","",ngs.name)
                fn = file.path(PGX.DIR,paste0(ngs.name1,".pgx"))
                dbg("upload_compute :: saving PGX to",fn)
                save(ngs, file=fn)  ## would clash with multiple users...
                ##pgx.dir=PGX.DIR;inc.progress=FALSE;pgx=PGXINFO
                pgx.name = paste0(sub("[.]pgx","",ngs$name),".pgx")        
                old.info = PGXINFO[which(!PGXINFO$dataset %in% c(ngs.name,pgx.name)),,drop=FALSE]
                new.info <- pgx.scanInfo(pgx.dir=PGX.DIR, pgx=old.info, verbose=FALSE)
                PGXINFO <<- new.info
                Sys.chmod(PGXINFO.FILE, mode="0666")
                write.csv(PGXINFO, file=PGXINFO.FILE)
            }
            showModal( modalDialog(
                "Saving data...",
                title = NULL, size = "s", fade=FALSE
            ))
            withProgress(message='saving data...', value=0.9, Sys.sleep(3))
            removeModal()
            
        } else {
            showModal( modalDialog(
                "Sorry. Saving only possible in Pro mode",
                title = NULL, size = "s", fade=FALSE
            ))
        }

    })


    uploadStatusTable <- reactive({
        
        dbg("[uploadStatusTable] uploaded_files$name=",input$upload_files$name)
        dbg("[uploadStatusTable] uploaded_files$datapath=",input$upload_files$datapath)

        ##for(i in 1:length(uploaded_files)) uploaded_files[[i]] <- NULL
        uploaded_files[["uploaded.pgx"]] <- NULL
        
        ## read uploaded files
        from.pgx = FALSE
        has.pgx <- any(grepl("[.]pgx$",input$upload_files$name))
        ff <- list()
        if(has.pgx) {
            i <- grep("[.]pgx$",input$upload_files$name)
            load(input$upload_files$datapath[i])  ## load NGS/PGX
            ff[["counts.csv"]] <- ngs$counts
            ff[["samples.csv"]] <- ngs$samples
            ff[["contrasts.csv"]] <- ngs$model.parameters$contr.matrix
            
            if(is.null(ngs$name)) ngs$name <- sub(".pgx$","",input$upload_files$name[i])
            uploaded_files[["uploaded.pgx"]] <- ngs
            from.pgx = TRUE

        } else {
            ii <- grep("csv$",input$upload_files$name)
            inputnames <- input$upload_files$name[ii]
            uploadnames <- input$upload_files$datapath[ii]
            dbg("[uploadStatusTable] inputnames=",inputnames,"\n")
            dbg("[uploadStatusTable] uploadnames=",uploadnames,"\n")
            if(length(uploadnames)>0) {
                for(i in 1:length(uploadnames)) {
                    fn1 <- inputnames[i]
                    fn2 <- uploadnames[i]
                    df <- NULL
                    if(grepl("counts",fn1)) {
                        ## allows duplicated rownames
                        df0 <- read.csv(fn2, check.names=FALSE, stringsAsFactors=FALSE)
                        df <- as.matrix(df0[,-1])
                        rownames(df) <- as.character(df0[,1])
                    } else {
                        df <- read.csv(fn2, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
                    }
                    ff[[ inputnames[i] ]] <- df
                }
            }            
        }

        
        ## store files in reactive value
        files.needed = c("counts.csv","samples.csv","contrasts.csv")
        ff = ff[ which(names(ff) %in% files.needed) ]
        if(length(ff)>0) {
            for(i in 1:length(ff)) {
                colnames(ff[[i]]) <- gsub("[\n\t ]","_",colnames(ff[[i]]))
                rownames(ff[[i]]) <- gsub("[\n\t ]","_",rownames(ff[[i]]))
                if(names(ff)[i] %in% c("counts.csv","contrasts.csv")) {
                    ff[[i]] <- as.matrix(ff[[i]])
                }
                uploaded_files[[names(ff)[i]]] <- ff[[i]]
            }
        }
        
        ## check dimensions
        files.uploaded <- names(uploaded_files)
        ##files.uploaded <- files.uploaded[match(files.needed,names(files.uploaded))]
        status = rep("please upload",3)
        names(status) = files.needed
        files.nrow = rep(NA,3)
        files.ncol = rep(NA,3)
        for(i in 1:3) {
            fn = files.needed[i]
            if(fn %in% files.uploaded) {
                status[i] = "OK"
                files.nrow[i] = nrow(uploaded_files[[fn]])
                files.ncol[i] = ncol(uploaded_files[[fn]])
            }
        }
        
        if(from.pgx==FALSE) {

            ## check files: matching dimensions
            if(status["counts.csv"]=="OK" && status["samples.csv"]=="OK") {
                if(!all( sort(colnames(uploaded_files[["counts.csv"]])) ==
                         sort(rownames(uploaded_files[["samples.csv"]])) )) {
                    status["counts.csv"] = "ERROR: colnames do not match (with samples)"
                    status["samples.csv"]  = "ERROR: rownames do not match (with counts)"
                }
            }
            
            MAXSAMPLES   = 25
            MAXCONTRASTS = 5
            MAXSAMPLES   = as.integer(max.limits["samples"])
            MAXCONTRASTS = as.integer(max.limits["comparisons"])
            
            ## check files: maximum contrasts allowed
            if(status["contrasts.csv"]=="OK") {
                if( ncol(uploaded_files[["contrasts.csv"]]) > MAXCONTRASTS ) {
                    status["contrasts.csv"] = paste("ERROR: max",MAXCONTRASTS,"contrasts allowed")
                }
            }
            
            ## check files: maximum samples allowed
            if(status["counts.csv"]=="OK" && status["samples.csv"]=="OK") {
                if( ncol(uploaded_files[["counts.csv"]]) > MAXSAMPLES ) {
                    status["counts.csv"]  = paste("ERROR: max",MAXSAMPLES," samples allowed")
                }
                if( nrow(uploaded_files[["samples.csv"]]) > MAXSAMPLES ) {
                    status["samples.csv"] = paste("ERROR: max",MAXSAMPLES,"samples allowed")
                }
            }
            
            ## check files: must have group column defined
            if(status["samples.csv"]=="OK" && status["contrasts.csv"]=="OK") {
                samples1   = uploaded_files[["samples.csv"]]
                contrasts1 = uploaded_files[["contrasts.csv"]]
                has.group <- "group" %in% colnames(samples1)
                matching.group <- all(rownames(contrasts1) %in% samples1$group)
                has.group
                matching.group
                cat("<uploaded_files> 3c : has.group=",has.group,"\n")
                cat("<uploaded_files> 3c : matching.group=",matching.group,"\n")
                if(!has.group) {
                    status["samples.csv"] = "ERROR: missing 'group' column"
                }
                if(has.group && !matching.group) {
                    status["contrasts.csv"] = "ERROR: contrasts do not match groups"
                }
            }
        } ## end-if-from-pgx
        
        ## check files
        description = c(
            "Count/expression file with gene on rows, samples as columns.",
            "Samples file with samples on rows, phenotypes as columns.",
            ## "Gene information file with genes on rows, gene info as columns.",
            "Contrast file with conditions on rows, contrasts as columns."        
        )
        df <- data.frame( status=status, filename=files.needed,
                         description=description,
                         nrow=files.nrow, ncol=files.ncol )

        ## deselect
        ## selectRows(proxy = dataTableProxy("pgxtable"), selected=NULL)
        return(df)    
    })

    output$uploadStatusTableOutput <- renderTable({
        uploadStatusTable()
    })

    userDatasetsTable <- reactive({
        df <- PGXINFO[1:2,] ## test-example
        ## df$description <- NULL
        df$datatype <- NULL
        df$conditions <- NULL        
        df$organism <- NULL
        df$nsets <- NULL
        df$path <- NULL
        ##df$sharing <- "private/public"
        
        df <- rbind(df,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)        
        df <- head(df,5)

        buttonInput <- function(FUN, len, id, ...) {
            inputs <- character(len)
            for (i in seq_len(len)) {
                inputs[i] <- as.character(FUN(paste0(id, i), ...))
            }
            ##inputs = paste("<div style='vertical-align: bottom;'",inputs,"</div>")
            ##inputs = paste("<div style='padding-top: 8px;'",inputs,"</div>")
            inputs
        }    
        sharing <- buttonInput(
            ##FUN = shinyWidgets::materialSwitch,
            ##FUN = shinyWidgets::switchInput,
            FUN = actionButton,
            len = nrow(df),
            id = 'uploaded_share_button_',
            ##size = "mini",
            width="80px",
            inline=TRUE,
            label = "publish",
            icon = icon("universal-access"),
            style='padding:2px; font-size:90%; color: black;'
        )
        delete <- buttonInput(
            FUN = actionButton,
            len = nrow(df),
            id = 'uploaded_delete_button_',
            label = "",
            ##size = "mini",
            width="50px",
            inline=TRUE,
            icon = icon("trash"),
            style='padding:2px; font-size:90%; color: #B22222;'
        )
        df <- cbind(sharing, delete, df)
        return(df)    
    })

    output$userDatasetsUI <- DT::renderDataTable({
        df <- userDatasetsTable()
        narrow.cols <- match(c("sharing","delete","nsamples","ngenes"),colnames(df))-1
        narrow.cols <- setdiff(narrow.cols,NA)
        DT::datatable(
                df,
                rownames=FALSE, escape = c(-1,-2),
                extensions = c('Scroller'),
                selection = list(mode='single', target='row', selected=1),                    
                fillContainer = TRUE,
                options = list(
                    ## dom = 'T<"clear">lfrtip',
                    dom = 'frti',
                    autoWidth = TRUE, ## scrollX=TRUE,
                    ordering = FALSE, paging=FALSE, searching=FALSE, info=FALSE,
                    scrollY = '100vh', ## scroller=TRUE,
                    columnDefs = list(list(width='5%', targets=narrow.cols))
                )
            ) %>%
            DT::formatStyle(0, target='row', fontSize='12px', lineHeight='90%')
    })


    ##---------------------------------------------------------------
    ##------------- modules for UsersMap ---------------------------
    ##---------------------------------------------------------------

    getUsersMapTable <- reactive({
        require(rgeolocate)

        access.files <- c(file.path(FILESX,"ncov2019_access.log"),
                          "/var/www/html/logs/access.log",
                          "/var/log/apache2/access.log")
        access.files <- access.files[file.exists(access.files)]
        access.files
        if(length(access.files)==0) return(NULL)
        ##accessfile = file.path(FILESX,"access-ncov2019.log")
        accessfile <- access.files[1]
        ##if(!file.exists(accessfile)) return(NULL)

        ## extract IP
        accessfile                          
        acc <- read.table(accessfile)
        ip <- as.character(acc[,1])
        ##loc <- ip_api(unique(ip))
        ip <- unique(ip)

        ## extract period
        acc.date <- gsub("[:].*|\\[","",as.character(acc[,4]))
        from.date <- head(acc.date,1)
        to.date <- tail(acc.date,1)
        from.to <- paste(from.date,"-",to.date)
        from.to
        
        ##file <- system.file("extdata","GeoLite2-Country.mmdb", package = "rgeolocate")
        ##loc <- maxmind(ip, file, "country_code")
        file <- file.path(FILESX,"GeoLite2-City.mmdb")
        loc <- rgeolocate::maxmind(ip, file, c("country_code", "country_name", "city_name"))
        country_code <- unique(loc$country_code)
        names(country_code) <- loc[match(country_code,loc$country_code),"country_name"]
        tt <- table(loc$country_name)
        df <- data.frame( country_name = names(tt),
                         country_code = country_code[names(tt)],
                         visitors = (as.integer(tt)))

        res <- list(table=df, period=from.to)

    })
    
    usersmap.RENDER %<a-% reactive({

        require(rworldmap)
        require(RColorBrewer)

        df <- ACCESS.LOG$table
        ##df <- getUsersMapTable()$table
        
        ##sPDF <- getMap()  
        ##mapCountryData(sPDF, nameColumnToPlot='continent')

        sPDF <- joinCountryData2Map(
            df,
            joinCode = "ISO2",
            nameJoinColumn = "country_code")
        
        par(mai=c(0,0.4,0.2,1),xaxs="i",yaxs="i")
        mapParams <- mapCountryData(
            sPDF, nameColumnToPlot="visitors",
            ##mapTitle = "Number of unique IPs",
            mapTitle = "", addLegend='FALSE',
            colourPalette = RColorBrewer::brewer.pal(9,"Blues"),
            numCats=9, catMethod="logFixedWidth")   
                   
        ##add a modified legend using the same initial parameters as mapCountryData
        do.call( addMapLegend,
                c(mapParams, labelFontSize = 0.85, legendWidth = 1.2, legendShrink = 0.5,
                  legendMar = 4, horizontal = FALSE, legendArgs = NULL, tcl = -0.5,
                  sigFigs = 4, digits = 3)
                )
        
    })
    
    usersmap_info = "<strong>Visitors map.</strong> The world map shows the number of users visiting this site by unique IP."
    
    callModule(
        plotModule,
        id = "usersmap", ## label="a", 
        plotlib = "baseplot",
        func = usersmap.RENDER,
        func2 = usersmap.RENDER, 
        info.text = usersmap_info,
        ##options = usersmap_options,
        pdf.width=12, pdf.height=7, pdf.pointsize=13,
        height = c(450,600), width = c('auto',1000), res=72,
        ##datacsv = enrich_getWordFreq,
        title = "Number of visitors by country"
    )

    ##usersmap_caption = "<b>(a)</b> <b>Geo locate.</b>"
    output$usersmapInfo <- renderUI({
        ##u <- getUsersMapTable()
        u <- ACCESS.LOG
        df <- u$table
        rownames(df) <-  df$country_name
        tot.users <- sum(df$visitors)
        freq <- df$visitors
        names(freq) <- df$country_name
        top.countries <- head(sort(freq,dec=TRUE),10)
        top.countriesTT <- paste("<li>",names(top.countries),top.countries,collapse=" ")
        
        HTML(
            "<b>Total visitors:</b>",tot.users,"<br><br>",
            "<b>Top 10 countries:</b><br><ol>",top.countriesTT,"</ol><br>",
            "<b>Period:</b><br>",u$period,"<br><br>"
        )
    })
    
    output$usersmap_UI <- renderUI({
        fillCol(
            height = 600,
            fillRow(
                flex = c(1,4.5),
                wellPanel( uiOutput(ns("usersmapInfo"))),
                plotWidget(ns("usersmap"))
            )
        )
    })


    ##---------------------------------------------------------------
    ##----------------- modules for Forum ---------------------------
    ##---------------------------------------------------------------
    
    output$forum <- renderUI({
        
        parenturl <- paste0(session$clientData$url_protocol,
                            "//",session$clientData$url_hostname,
                            ":",session$clientData$url_port,
                            session$clientData$url_pathname)
        ## parenturl <- gsub("localhost","127.0.0.1",parenturl)
        parenturl <- URLencode(parenturl, TRUE)
        cat("[LoadingModule:forum] parenturl =",parenturl,"\n")
        src = paste0('https://groups.google.com/forum/embed/?place=forum/omicsplayground',
                     '&showsearch=true&showpopout=true&parenturl=',parenturl)
        cat("src = ",src,"\n")
        tags$iframe(id="forum_embed", src=src, height=600, width='100%',
                    ##seamless="seamless",
                    frameborder='no')
        ##HTML(src)
    })
         
    output$tweet <- renderUI({
        ## NOT WORKING YET...
        tags$a(class="twitter-timeline",
               href="https://twitter.com/bigomics?ref_src=twsrc%5Etfw")
        ##tags$script('twttr.widgets.load(document.getElementById("tweet"));')
    })
            
    output$forum_UI <- renderUI({
        fillCol(
            height = 550,
            fillRow(
                flex=c(4,0),
                htmlOutput(ns("forum"))
                ##uiOutput("tweet")
            )
        )
    })
    
    ##------------------------------------------------
    ## Module return object
    ##------------------------------------------------
    res <- list(
        inputData = inputData,
        ##inputData = currentPGX,
        usermode = reactive({ USERMODE() })
    )
    return(res)
}


