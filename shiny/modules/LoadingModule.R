
LoadingInputs <- function(id) {
    ns <- NS(id)  ## namespace
    tagList(
        uiOutput(ns("description")),
        uiOutput(ns("inputsUI"))
        ##uiOutput(ns("socialButtons"))
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
            tabPanel("Public datasets",uiOutput(ns("pgxtable_UI"))),
            tabPanel("Upload data",uiOutput(ns("upload_UI")))
        )
    )
}


LoadingModule <- function(input, output, session, hideUserMode=FALSE)
{
    ns <- session$ns ## NAMESPACE
    
    ##useShinyjs(rmd=TRUE)
    useShinyjs()
    useSweetAlert()
    SHOWSPLASH=TRUE
    ## SHOWSPLASH=FALSE
    
    LOGIN_AUTHENTICATION = "none"
    ##LOGIN_AUTHENTICATION = "register"
    ##LOGIN_AUTHENTICATION = "password"
    CREDENTIALS = NULL
    if(file.exists("CREDENTIALS")) {
        CREDENTIALS <- read.table("CREDENTIALS",row.names=1,
                                  header=TRUE, stringsAsFactors=FALSE)
        LOGIN_AUTHENTICATION = "password"
    }

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
    USERMODE <- reactiveVal( factor("BASIC",levels=USERLEVELS) )

    output$inputsUI <- renderUI({

        usermodeUI <- tipify(radioGroupButtons(
            inputId = ns("main_usermode"),
            label = "User mode:",
            choices = USERLEVELS,
            selected = "BASIC",
            status = "warning",
            checkIcon = list(yes = icon("ok", lib = "glyphicon"))),
            usermode_infotip, placement="bottom", options = list(container="body")
            )
        
        if(hideUserMode) usermodeUI <- NULL
        
        ui <- tagList(
            usermodeUI,
            br(),br(),br(),
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

    observeEvent( input$action_beer, {
        dbg("buy beer button action\n")
        startup_count <<- startup_count + 1    
        USER$logged <- TRUE
        USER$name   <- "beer buddy"
        removeModal()
        ##alert("Wow. Thanks buddy!")
        sendSweetAlert(
            session=session, title="Wow. Thanks buddy!",
            text = "Free entrance for you!", type = "info")
        ##Sys.sleep(4);removeModal()
    })


    observeEvent( input$action_play, {

        dbg("action_play:: play button action\n")
        cat("action_play:: LOGIN_AUTHENTICATION=",LOGIN_AUTHENTICATION,"\n")

        startup_count <<- startup_count + 1    
        if(LOGIN_AUTHENTICATION!="none") {
            showLogin()  ## $
        } else {
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

    USER <- reactiveValues( logged = FALSE, name="anonymous")

    observeEvent( input$logout, {
        ##updateTextInput(session, ".username", value=NULL)
        reset(ns("login_username"))
        reset(ns("login_password"))
        USER$logged <- FALSE
    })

    if(LOGIN_AUTHENTICATION=="password") {

        showLogin <- function() {
            showModal( modalDialog(
                id = "login",
                title = "Login to Omics Playground",
                tagList( 
                    textInput(ns("login_username"), "Username:"),
                    passwordInput(ns("login_password"), "Password:"),
                    div( actionButton(ns("login_btn"), "Login"),
                        style="text-align: center;")
                ),
                footer = textOutput(ns("login_warning")),
                size = "s"
            ))
        }


    } else if(LOGIN_AUTHENTICATION=="register") {    

        showLogin <- function() {
            showModal( modalDialog(
                id = "login",
                title = "Login to Omics Playground",
                tagList( 
                    p("Please", actionLink(ns("register_link"),"register"),"to use the Omics Playground. If you have already registered please enter your registration email address below."),
                    textInput(ns("login_username"), "E-mail:"),
                    div( actionButton(ns("login_btn"), "Login"),
                        style="text-align: center;")
                ),
                footer = div(textOutput(ns("login_warning")),style="color: red;"),
                size = "s"
            ))
        }
        
        observeEvent( input$register_link, {
            ##install.packages("countrycode")
            require(countrycode)
            all.countries = countrycode::codelist[,"country.name.en"]
            
            showModal( modalDialog(
                id = "register",
                title = "Omics Playground registration",
                tagList( 
                    p("Fill in the form below. Registration data is used to measure usage only. This helps us track and better serve our user community."),
                    ## textInput("register_name", "Name:"),
                    textInput("register_email", "E-mail:"),
                    ## textInput("register_organization", "Organization:"),
                    selectInput("register_country", "Country:", choices=c("",all.countries)),
                    div( ##actionButton(ns("register_btn_skip"), "Skip"),
                        actionButton(ns("register_btn"), "Register"),
                        style="text-align: center;")
                ),
                footer = div(textOutput("register_warning"),style="color: red;"),
                size = "s"
            ))
        })           

        observeEvent( input$register_btn, {
            register.OK = ( input$register_email != "" &&
                            grepl("@",input$register_email) &&  ## valid email
                            ## input$register_name != "" &&
                            input$register_country != "" )
            if(register.OK) {
                rdata <- paste("Tom Cruise","ACME Inc.","USA",date(),sep=",")
                rdata <- paste(input$register_email,
                               input$register_name,
                               input$register_organization,
                               input$register_country,
                               date(), sep=",")
                write( rdata, file="logs/registered.csv", append=TRUE )

                USER$name   <- input$register_email
                USER$logged <- TRUE  ## global reactive
                show("register_warning")
                output$register_warning = renderText(paste("Welcome",input$register_email,"!"))
                delay(3000, hide("register_warning", anim = TRUE, animType = "fade"))
                removeModal()
                
            } else {
                ##show("register_warning")
                output$register_warning = renderText("Invalid email or country")
                ##delay(2000, hide("register_warning", anim = TRUE, animType = "fade"))
                delay(2000, {output$register_warning <- renderText("")})
            }
        })
    } else {
        showLogin <- function() {}
    }

    output$login_warning = renderText("")

    observeEvent( input$login_btn, {           

        cat("LOGIN_AUTHENTICATION=",LOGIN_AUTHENTICATION,"\n")
        cat("username=",input$login_username,"\n")
        cat("password=",input$login_password,"\n")
        cat("Logged=",USER$logged,"\n")

        login.OK = FALSE
        if(LOGIN_AUTHENTICATION=="register") {
            ##if( is.null(input$login_username) || is.null(input$login_password)) return(NULL)
            ##if( input$login_username=="" || input$login_password=="") return(NULL)    
            ##if( is.null(input$login_name) || input$login_name == "") return(NULL)
            if( is.null(input$login_username) || input$login_username=="") return(NULL)
            registered <- read.csv(file="logs/registered.csv",header=FALSE,stringsAsFactors=FALSE)
            registered.users <- unique(registered[,1]) ## emails

            cat("registered.users=",registered.users,"\n")
            registered.users <- c(registered.users,"demo@bigomics.ch")
            login.OK = (input$login_username %in% registered.users)
        }
        
        if(LOGIN_AUTHENTICATION=="password") {
            if( is.null(input$login_username) || is.null(input$login_password)) return(NULL)
            if( input$login_username=="" || input$login_password=="") return(NULL)    
            username <- input$login_username
            ok.user <- isTRUE(CREDENTIALS[username,"password"]==input$login_password)
            ok.date <- isTRUE( Sys.Date() < as.Date(CREDENTIALS[username,"expiry"]) )
            login.OK = (ok.user && ok.date)
        }

        if (login.OK) {
            output$login_warning = renderText("")
            removeModal()

            USER$name   <- input$login_username
            USER$logged <- TRUE
            
            ## Here you can perform some user-specific functions, or site news
            if(0) {
                showModal(modalDialog(
                    paste("Welcome",USER$name,"!"),
                    footer=NULL, size="m"))
                Sys.sleep(3)
            }
            removeModal()
            
        } else {
            ##show("login_warning")
            if(LOGIN_AUTHENTICATION=="password") {
                output$login_warning = renderText("Invalid username or password")
            }
            if(LOGIN_AUTHENTICATION=="register") {
                output$login_warning = renderText("Email address not recognized")
            }
            ##delay(2000, hide("login_warning", anim = TRUE, animType = "fade"))
            delay(2000, {output$login_warning <- renderText("")})
            USER$logged <- FALSE
        }
        ##hide("login_warning")
    })
    
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
        dbg("inputData:: ---------- reacted ---------------\n")
        dbg("inputData:: LOGIN_AUTHENTICATION=",LOGIN_AUTHENTICATION,"\n")
        dbg("inputData:: USER$Logged=",USER$logged,"\n")
        
        ## authenicate user
        ##if(!USER$logged) showLogin()
        if(LOGIN_AUTHENTICATION!="none" && !USER$logged) return(NULL)

        pgx <- currentPGX()

        return(pgx)
    })

    observeEvent( input$loadbutton, {

        ## Observe button press
        btn <- isolate(input$loadbutton)
        pgx = NULL
        pgx = isolate(selectedDataSet())

        dbg("[observe:loadbutton] loadbutton=",btn,"\n")
        dbg("[observe:loadbutton] 1: pgx.selected=",pgx)

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
            dbg("[observe:loadbutton] LOADING",pgx1,"\n")
            ##withProgress(message='loading...', value=0.8,
            load(pgx1,verbose=0)
        } else {
            cat("[observe:loadbutton] ERROR file not found : ",pgx1,"\n")
            removeModal()
            return(NULL)
        }
        
        ##----------------- update input
        dbg("[observe:loadbutton] head.names.PGX=",head(names(pgx)))
        dbg("[observe:loadbutton] initializing PGX object")
        ngs <- pgx.initialize(ngs)
        if(is.null(ngs)) {
            cat("[observe:loadbutton] ERROR in object initialization\n")
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
        dbg("[observe:loadbutton] ready! \n")
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

    ## Some 'global' reactive variables used in this file
    uploaded_files <- reactiveValues()
    ##shinyjs::disable("upload_compute")

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
        
        dbg("<getPGXTable> 2")
        
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
            flex=c(1),
            height = 580,
            ##moduleWidget(pgxplots_module, outputFunc="plotOutput")
            ##moduleWidget(pgxtable_module, outputFunc="dataTableOutput", ns=ns)
            tableWidget(ns("pgxtable"))
        )
    })

    ##================================================================================
    ## Upload data
    ##================================================================================
    
    output$downloadExampleData <- downloadHandler(
        filename = "exampledata.zip",
        content = function(file) {
            zip = file.path(FILES,"exampledata.zip")
            file.copy(zip,file)
        }
    )

    upload_info = "<h4>User file upload</h4><p>Please prepare the data files in CSV format as listed below. It is important to name the files exactly as shown. The file format must be comma-separated-values (CSV) text. Be sure the dimensions, rownames and column names match for all files. You can download a zip file with example files here: EXAMPLEZIP. You can upload a maximum of <u>LIMITS</u>. If you want to analyze larger datasets, please use the scripts and upload the pgx file. After uploading, in sidebar on the left, provide a name for your dataset. Finally, hit the compute button. The computations may take 10 to 30 minutes depending on the size of your dataset."
    DLlink = downloadLink(ns("downloadExampleData"),"exampledata.zip")
    upload_info = sub("EXAMPLEZIP", DLlink, upload_info)
    
    upload_info2 =
        "<br><h4>Uploaded datasets</h4><p>Below are your uploaded datasets. As a free user, you can only have a maximum of one private dataset. If you want to analyze a new dataset, you must either delete your old dataset or make the dataset public.<br><br>"

    ##upload_filetypes = c("text/csv","text/comma-separated-values,text/plain",".csv")
    upload_filetypes = c(".csv",".pgx")    
    output$upload_UI <- renderUI({    

        basic.limits = "25 samples and 5 comparisons"
        pro.limits   = "1000 samples and 20 comparisons"
        if(USERMODE()=="BASIC") upload_info = sub("LIMITS", basic.limits, upload_info)
        if(USERMODE()=="PRO") upload_info = sub("LIMITS", pro.limits, upload_info)    
        userdataUI <- NULL
        if(DEV.VERSION) {
            userdataUI <- DT::dataTableOutput(ns("userDatasetsUI"))
        }
        
        fillCol(
            flex = c(1,1.5),
            height = 750,
            fillRow(
                flex = c(1,0.1,3.5),
                wellPanel(
                    fileInput(ns("upload_files"), "Choose files",
                              multiple = TRUE, accept = upload_filetypes),
                    textInput(ns("upload_name"),"Name of dataset:"),
                    ##textAreaInput("upload_description", "Description:", value = NULL,
                    ##              rows=5, placeholder="Describe your data set (minimum 100 characters)"),
                    actionButton(ns("upload_compute"),"Compute!",icon=icon("running"))
                ),br(),
                fillCol(
                    flex = c(NA,1),
                    div(HTML(upload_info),style="font-size: 13px;"),
                    tableOutput(ns("upload_status"))
                )
            ),
            fillRow(
                userdataUI
            )
        )
    })

    output$allFilesOK <- reactive({    
        files.needed = c("counts.csv","samples.csv","contrasts.csv")
        all.there <- all(files.needed %in% names(uploaded_files))    
        df <- uploadStatusTable()
        filled = all.ok = TRUE
        all.ok = all( df$status == "OK")
        filled = (input$upload_name!="")
        ## filled = (input$upload_name!="" && nchar(input$upload_description)>=100)
        ## all.ok    <- all(sapply(uploaded_files, function(x) x$status=="OK"))
        ok <- (all.there && all.ok && filled)
        if(ok) shinyjs::enable("upload_compute")
        if(!ok) shinyjs::disable("upload_compute")
        ok
    })
    outputOptions(output, "allFilesOK", suspendWhenHidden = FALSE) ## important!

    observeEvent( input$upload_compute, {
        ## are you sure? Any message/warning before computation is done.
        if(0) {
            require(shinyWidgets)
            confirmSweetAlert(
                session = session,
                inputId = ns("myconfirmation"),
                text = "Your data will be uploaded for computation. Are you sure?"
            )
        } else {
            showModal( modalDialog(
                HTML("Your data will be uploaded for computation. By uploading the data you accept our EULA and you confirm that the data has been anonymized. For research use only."),
                title = NULL,
                size = "s",
                footer = tagList(
                    modalButton("Cancel"),
                    actionButton(ns("myconfirmation2"),"Confirm", icon=NULL)
                )
            ))
        }
    })

    observeEvent( c(input$myconfirmation,input$myconfirmation2), {

        ## if(isFALSE(input$myconfirmation)) { return(NULL) }

        ## --------------------- OK start ---------------------------
        dbg("upload_compute :: showing coffee modal")
        showLoadingModal("Calculating... it's a good time to get a coffee now.")

        has.pgx <- ("uploaded.pgx" %in% names(uploaded_files))       
        has.pgx <- has.pgx && !is.null(uploaded_files[["uploaded.pgx"]])
        cat("upload_compute: names(uploaded_files)=",names(uploaded_files),"\n")
        
        if(has.pgx) {

            dbg("upload_compute :: ***** using 'uploaded.pgx' ******")        
            ngs <- uploaded_files[["uploaded.pgx"]]
            
        } else {
            dbg("upload_compute :: ***** real computation *****")
            
            counts    <- as.matrix(uploaded_files[["counts.csv"]])
            samples   <- data.frame(uploaded_files[["samples.csv"]],stringsAsFactors=FALSE)
            contrasts <- as.matrix(uploaded_files[["contrasts.csv"]])

            max.genes = NULL
            if( USERMODE() == "BASIC") {
                gx.methods   = c("ttest.welch","ttest.rank","trend.limma") ## fastest 3
                gset.methods = c("fisher","gsva","camera")  ## fastest 3            
                ## gx.methods   = c("trend.limma","edger.qlf","edger.lrt")
                ## gset.methods = c("fisher","gsva","fgsea")
                extra.methods = c("meta.go","infer","drugs","wordcloud")
                max.genes = 10000
            } else {
                gx.methods   = c("ttest.welch","trend.limma","edger.qlf","deseq2.wald")
                gset.methods = c("fisher","gsva","fgsea","camera","fry")
                extra.methods = c("meta.go","infer","deconv","drugs-combo","wordcloud")
                max.genes = 25000
                if(ncol(counts) > 1000) {
                    ## probably scRNA-seq
                    gx.methods   = c("ttest","ttest.welch","trend.limma")
                    gset.methods = c("fisher","gsva","fgsea")
                    extra.methods = c("meta.go","infer","deconv","drugs-combo","wordcloud")                    
                    max.genes = 10000
                }
            }
            if(DEV.VERSION) {
                gx.methods   = c("ttest","ttest.rank","ttest.welch","trend.limma","edger.qlf","edger.lrt","deseq2.wald")
                gset.methods = c("fisher","gsva","fgsea","camera","fry","ssgsea","spearman")
                extra.methods = c("meta.go","infer","deconv","drugs-combo","wordcloud")
                max.genes = 9999999
            }
            
            ##extra.methods = c("meta.go","infer")
            
            ##----------------------------------------------------------------------
            ## Upload and do precomputation
            ##----------------------------------------------------------------------
            start_time <- Sys.time()
            ## Create a Progress object
            progress <- shiny::Progress$new()
            on.exit(progress$close())    
            progress$set(message = "Processing", value = 0)

            ngs <- pgx.upload(
                counts, samples, contrasts,
                max.genes = max.genes,
                progress = progress,
                gx.methods = gx.methods,
                gset.methods = gset.methods,
                extra.methods = extra.methods,
                lib.dir = FILES,
                only.hugo = TRUE
            )

            end_time <- Sys.time()
            delta_time  = end_time - start_time
            delta_time
            dbg("upload_compute :: total processing time of",delta_time,"secs")
            
            names(ngs)
            head(ngs$samples)
            
            ngs$name = "uploaded"
            ngs$datatype = "generic"
            ngs$description = "not available"
            
            ngs.name = gsub("[ ]","-",input$upload_name)
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
                modalButton("Dismiss")
            )
        ))
        
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
        dbg("[home] observeEvent:sharedata\n")
        showModal( modalDialog(
            HTML("<center><b>New feature!</b><br>Sharing your data with others will be soon available as new feature!</center>"),
            title = NULL, size = "s", fade=FALSE
        ))
        return(NULL)
    })
    
    observeEvent( input$savedata, {
        dbg("[home] observeEvent:savedata\n")

        showModal( modalDialog(
            HTML("<center><b>New feature!</b><br>Saving your data to the cloud will be soon available as new feature for Pro users!</center>"),
            title = NULL, size = "s", fade=FALSE
        ))
        return(NULL)
        
        ## -------------------- save PGX file/object
        saving.ok = ( USERMODE() == "PRO")
        if(saving.ok) {

            if(0) {
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
        
        cat("<uploaded_files> name=",input$upload_files$name,"\n")
        cat("<uploaded_files> datapath=",input$upload_files$datapath,"\n")
        ##for(i in 1:length(uploaded_files)) uploaded_files[[i]] <- NULL
        uploaded_files[["uploaded.pgx"]] <- NULL
        
        ## read uploaded files
        from.pgx = FALSE
        has.pgx <- any(grepl("[.]pgx$",input$upload_files$name))
        if(has.pgx) {
            i <- grep("[.]pgx$",input$upload_files$name)
            load(input$upload_files$datapath[i])  ## load NGS/PGX
            ff <- list()
            ff[["counts.csv"]] <- ngs$counts
            ff[["samples.csv"]] <- ngs$samples
            ff[["contrasts.csv"]] <- ngs$model.parameters$contr.matrix

            if(is.null(ngs$name)) ngs$name <- sub(".pgx$","",input$upload_files$name[i])
            uploaded_files[["uploaded.pgx"]] <- ngs
            from.pgx = TRUE
        } else {
            ii <- grep("csv$",input$upload_files$name)
            ff = lapply(input$upload_files$datapath[ii], read.csv, row.names=1,
                        check.names=FALSE, stringsAsFactors=FALSE )
            names(ff) <- input$upload_files$name[ii]    
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

        if(!from.pgx) {

            ## check files: matching dimensions
            if(status["counts.csv"]=="OK" && status["samples.csv"]=="OK") {
                if(!all( colnames(uploaded_files[["counts.csv"]]) ==
                         rownames(uploaded_files[["samples.csv"]]) )) {
                    status["counts.csv"] = "ERROR: colnames do not match (with samples)"
                    status["samples.csv"]  = "ERROR: rownames do not match (with counts)"
                }
            }
            
            MAXSAMPLES   = 25
            MAXCONTRASTS = 5

            if( USERMODE() == "PRO") {
                MAXSAMPLES   = 1000
                MAXCONTRASTS = 20
            }
            if( USERMODE() == "DEV") {
                MAXSAMPLES   = 999999
                MAXCONTRASTS = 99999
            }
            
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
        }

        
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

    output$upload_status <- renderTable({
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


    ##------------------------------------------------
    ## Module return object
    ##------------------------------------------------
    res <- list(
        inputData = inputData,
        usermode = reactive({ USERMODE() })
    )
    return(res)
}


