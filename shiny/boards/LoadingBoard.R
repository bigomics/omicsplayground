##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing LoadingBoard")

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
        tabsetPanel(
            id = ns("tabs"),
            tabPanel("Datasets",uiOutput(ns("pgxtable_UI"))),
            tabPanel("Upload data",uiOutput(ns("upload_UI"))),
            tabPanel("Visitors map",uiOutput(ns("usersmap_UI")))
            ## tabPanel("Community forum",uiOutput(ns("forum_UI")))
        )
    )
}

LoadingBoard <- function(input, output, session, 
                         max.limits = c("samples"=1000,"comparisons"=20,
                                        "genes"=20000, "genesets"=10000),
                         enable_delete = TRUE, enable_save = TRUE,
                         authentication="none", firebase=NULL, firebase2=NULL)
{
    ns <- session$ns ## NAMESPACE
    
    ##useShinyjs(rmd=TRUE)
    useShinyjs()
    useSweetAlert()
    SHOWSPLASH=TRUE
    ## SHOWSPLASH=FALSE

    message("[LoadingBoard] in.shinyproxy = ",in.shinyproxy())    
    message("[LoadingBoard] SHINYPROXY_USERNAME = ",Sys.getenv("SHINYPROXY_USERNAME"))
    message("[LoadingBoard] SHINYPROXY_USERGROUPS = ",Sys.getenv("SHINYPROXY_USERGROUPS"))

    AUTHENTICATION = authentication
    message("[LoadingBoard] AUTHENTICATION=",AUTHENTICATION)

    auth <- NULL
    if(AUTHENTICATION == "password") {
        auth <- callModule(
            PasswordAuthenticationModule, "auth",
            credentials.file = "CREDENTIALS")
    } else if(AUTHENTICATION == "firebase") {
        require(firebase)
        ##firebase <- FirebaseEmailPassword$new()
        auth <- callModule(
            ##FirebaseAuthenticationModule, "auth")
            FirebaseAuthenticationModule, "auth",
            firebase = firebase, firebase2 = firebase2)
    } else if(AUTHENTICATION == "register") {
        auth <- callModule(
            RegisterAuthenticationModule, "auth",
            register.file = "../logs/register.log")
    } else if(AUTHENTICATION == "shinyproxy" && in.shinyproxy()) {
        username <- NULL
        is.anonymous <- Sys.getenv("SHINYPROXY_USERGROUPS")=="ANONYMOUS"
        if(!is.anonymous) username <- Sys.getenv("SHINYPROXY_USERNAME")
        auth <- callModule(NoAuthenticationModule, "auth", username=username)
    } else {
        ## none
        auth <- callModule(NoAuthenticationModule, "auth")
    } 
    
    ##-----------------------------------------------------------------------------
    ## Description
    ##-----------------------------------------------------------------------------
    description = "<b>Omics Playground</b> is a self-service bioinformatics platform for interactive analysis, visualization and interpretation of transcriptomics and proteomics data. Life scientists can easily perform complex data analysis and visualization without coding, and significantly reduce the time to discovery."
    
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
    downloadButton2 <- function (outputId, label = "Download", class = NULL, ...) {
        aTag <- tags$a(id = outputId,
                       class = paste("btn btn-default shiny-download-link", class),
                       href = "", target = "_blank", download = NA, 
                       icon("file-csv"), label, ...)
    }
    
    output$inputsUI <- renderUI({        

        delete_button <- NULL
        if(enable_delete) {
            delete_button <- tipify( actionButton(
                ns("deletebutton"), label=NULL, icon=icon("trash"),
                style='padding:2px 1px 1px 1px; font-size:140%; color: #B22222; width:30px;'
            ),"Delete the selected dataset.", placement="bottom")
        }

        ui <- tagList(
            useShinyalert(),  # Set up shinyalert
            p(strong("Dataset info:")),
            div( htmlOutput(ns("dataset_info")), id="datainfo"),
            br(),
            conditionalPanel(
                "output.rowselected != 0", ns=ns,
                tipify( actionButton(ns("loadbutton"),label="Load",class="load-button"),
                   "Click to load the selected dataset.", placement="bottom"),
                tipify( downloadButton(
                    ns("downloadpgx"), label=NULL, ## icon=icon("download"),
                    style='padding:2px 1px 1px 1px; font-size:140%; width:30px;'
                ),"Download PGX file (binary).", placement="bottom"),
                tipify( downloadButton2(
                    ns("downloadzip"), label=NULL, icon=icon("file-csv"),
                    style='padding:2px 1px 1px 1px; font-size:140%; width:30px;'
                ),"Download CSV files (counts.csv, samples.csv, contrasts.csv).",
                placement="bottom"),
                delete_button
            ),
            br(),br(),
            tipify( actionLink(ns("showfilter"), "show filters", icon=icon("cog", lib = "glyphicon")),
                   "Show dataset filters.", placement="top"),
            br(),br(),
            conditionalPanel(
                "input.showfilter % 2 == 1", ns=ns,
                uiOutput(ns("dataset_filter"))
            )
        )
        ui
    })
    outputOptions(output, "inputsUI", suspendWhenHidden=FALSE)
    output$rowselected <- reactive({ !is.null(selectedPGX()) })
    outputOptions(output, "rowselected", suspendWhenHidden=FALSE)
    
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
            HTML(paste("<p><b>",names(inf),"</b>:", inf,""))
        }
    })

    output$dataset_filter <- renderUI({
        df <- PGXINFO()
        collections <- sort(setdiff(df$collection,c(NA,"")))
        datatypes <- sort(setdiff(df$datatype,c(NA,"")))
        organisms <- sort(setdiff(df$organism,c(NA,"")))        
        tagList(
            checkboxGroupInput(ns("flt_datasets"),"datasets", choices = collections),
            checkboxGroupInput(ns("flt_datatype"),"datatype", choices = datatypes),
            checkboxGroupInput(ns("flt_organism"),"organism", choices = organisms)
            ##checkboxGroupInput("flt_conditions","conditions",
            ##choices=c("treatment","sex","activated"))
        )

    })
    
    ##-----------------------------------------------------------------------------
    ## READ initial PGX file info
    ##-----------------------------------------------------------------------------
    ##PGXINFO <- pgx.updateInfoFile(PGX.DIR, file="datasets-info.csv", 
    ##                           force=FALSE, verbose=TRUE )
    PGXINFO.FILE <- file.path(PGX.DIR[1], "datasets-info.csv")  ## first folder!!
    PGXINFO  <- reactiveVal(NULL)    
    infofile <- pgx.scanInfoFile(PGX.DIR, file="datasets-info.csv", verbose=TRUE )
    cat("[LoadingBoard] dim.infofile = ",dim(infofile),"\n")
    if(!is.null(infofile) && nrow(infofile)) {
        if(!"collection" %in% colnames(infofile)) {
            infofile$collection <- rep("",nrow(infofile))
        } else {
            infofile$collection <- ifelse(is.na(infofile$collection),"",infofile$collection)
        }
    }
    PGXINFO(infofile)

    selectedPGX <- reactive({
        ##sel <- input$pgxtable_rows_selected
        sel <- pgxtable$rows_selected()
        if(is.null(sel) || length(sel)==0) return(NULL)
        df <- getPGXTable()
        pgxfile <- as.character(df$dataset[sel])
        pgxfile <- paste0(sub("[.]pgx$","",pgxfile),".pgx") ## add/replace .pgx
        pgxfile
    })

    ##=============================================================================
    ##========================== OBSERVE/REACT ====================================
    ##=============================================================================
    pgxfile="geiger2016-arginine"

    loadPGX <- function(pgxfile) {        
        pgxfile <- paste0(sub("[.]pgx$","",pgxfile),".pgx") ## add/replace .pgx         
        pgx.path <- PGX.DIR[file.exists(file.path(PGX.DIR,pgxfile))][1]
        pgxfile1 = file.path(pgx.path,pgxfile)
        pgxfile1
        ngs <- NULL
        pgx <- NULL
        if(file.exists(pgxfile1)) {
            message("[LoadingBoard::loadPGX] loading ",pgxfile)
            ##withProgress(message='loading...', value=0.8,
            load(pgxfile1,verbose=0)
        } else {
            cat("[LoadingBoard::loadPGX] error file not found : ",pgxfile)
            return(NULL)
        }
        if(!is.null(ngs)) return(ngs)
        if(!is.null(pgx)) return(pgx)
    }
    
    output$downloadpgx <- downloadHandler(
        ##filename = "userdata.pgx",
        filename = function() {
            selectedPGX()
        },
        content = function(file) {
            pgxfile <- selectedPGX()
            cat("[LoadingBoard::loadPGX] pgxfile = ",pgxfile,"\n")
            if(is.null(pgxfile) || pgxfile=="" || length(pgxfile)==0) return(NULL)
            ngs <- loadPGX(pgxfile)
            temp <- tempfile()
            cat("[LoadingBoard::loadPGX] temp = ",temp)            
            save(ngs, file=temp)
            file.copy(temp,file)
        }
    )

    output$downloadzip <- downloadHandler(
        ##filename = "userdata.zip",
        filename = function() {
            sub("pgx$","zip",selectedPGX())
        },
        content = function(file) {
            ## ngs <- currentPGX()  ## current dataset
            pgxfile <- selectedPGX()
            cat("[LoadingBoard::downloadZIP] pgxfile = ",pgxfile,"\n")            
            if(is.null(pgxfile) || pgxfile=="" || length(pgxfile)==0) return(NULL)
            pgxname <- sub("[.]pgx$","",pgxfile)
            ngs <- loadPGX(pgxfile)
            dir.create(tmp <- tempfile())
            tmp2 <- file.path(tmp,pgxname)
            dir.create(tmp2)
            write.csv(ngs$counts,  file=file.path(tmp2, "counts.csv"))
            write.csv(ngs$samples, file=file.path(tmp2, "samples.csv"))
            write.csv(ngs$model.parameters$exp.matrix, file=file.path(tmp2, "contrasts.csv"))
            zipfile <- tempfile(fileext = ".zip")
            zip::zip(zipfile,
                     files=paste0(pgxname,"/",c("counts.csv","samples.csv","contrasts.csv")),
                     root=tmp)
            ## zip::zip_list(zipfile)
            cat("[LoadingBoard::downloadZIP] zipfile = ",zipfile)                        
            file.copy(zipfile,file)
        }
    )

    touchtable <- reactiveVal(0)
    
    observeEvent( input$deletebutton, {
                
        ##pgxfile <- currentPGX()$name
        pgxfile <- selectedPGX()
        if(is.null(pgxfile) || pgxfile=="" || length(pgxfile)==0) return(NULL)
        pgx.path <- PGX.DIR[file.exists(file.path(PGX.DIR,pgxfile))][1]
        pgxfile1 = file.path(pgx.path,pgxfile)
        pgxfile1
        sel <- NULL

        deletePGX <- function() {
            if(input$confirmdelete) {

                cat(">>> deleting",pgxfile,"\n")
                pgxfile2 <- paste0(pgxfile1,"_")  ## mark as deleted
                file.rename(pgxfile1, pgxfile2)
                ##touchtable(touchtable()+1)

                this.pgx <- sub("[.]pgx$","",pgxfile)
                all.pgx  <- sub("[.]pgx$","",PGXINFO()$dataset)

                ## get selected row before deleting
                table.pgx <- sub("[.]pgx$","",getPGXTable()$dataset)
                sel <- which(table.pgx == this.pgx)

                newpgx <- PGXINFO()[all.pgx != this.pgx,]
                PGXINFO(newpgx)

                selectRows(proxy=dataTableProxy(ns("pgxtable")), selected=sel)
            } else {
                cat(">>> deletion cancelled\n")
            }
        }


        sel  <- which(sub("[.]pgx$","",PGXINFO()$dataset) == sub("[.]pgx$","",pgxfile))
        this.pgxinfo <- PGXINFO()[sel,]

        owner1 = "owner"
        owner1 <- this.pgxinfo$owner
        if(is.null(owner1) || is.na(owner1)) owner1 <- ""
        is.owner <- (owner1 == auth$name())
        ## must be owner and not empty/anonymous
        not.anonymous <- !is.na(auth$name()) && auth$name()!="" 
        
        is.uploaded <- this.pgxinfo$collection == "uploaded"
        ##allow.delete <- is.owner && !is.na(owner1) && owner1!="" && 
        ##    !is.na(auth$name() && auth$name()!="" )
        allow.delete <- is.uploaded && !not.anonymous
        
        message("[LoadingBoard::@deletebutton] WARNING:: ",pgxfile," owned by ",owner1," \n")
        message("[LoadingBoard::@deletebutton] current user = ",auth$name()," \n")
        message("[LoadingBoard::@deletebutton] allow.delete = ",allow.delete," \n")
        message("[LoadingBoard::@deletebutton] is.owner = ",is.owner," \n")
        
        allow.delete = TRUE
        if(!allow.delete) {
            message("[LoadingBoard::@deletebutton] WARNING:: ",pgxfile," not owned by ",auth$name()," \n")
            shinyalert::shinyalert(
                            title = "Error!",
                            text = "You do not have permission to delete this dataset",
                            type = "error"
                        )
        } else {
            shinyalert::shinyalert(
                "Delete this dataset?",
                paste("Are you sure you want\nto delete '",pgxfile,"'?"),
                confirmButtonText = "Delete",
                showCancelButton = TRUE,
                callbackR = deletePGX,
                inputId = "confirmdelete")
        }

        
    })
    
    ##=================================================================================
    ##=============================== MODAL DIALOGS ===================================
    ##=================================================================================

    
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
        AuthenticationUI(ns("auth"))
        
        startup_count <<- startup_count + 1    
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

    selectedDataSetInfo <- reactive({
        ##sel <- input$pgxtable_rows_selected
        sel <- pgxtable$rows_selected()
        if(is.null(sel) || length(sel)==0) return(NULL)
        df <- getPGXTable()
        unlist(lapply(df[sel,],as.character))
    })

    currentSection <- reactive({
        cdata <- session$clientData
        sub("section-","",cdata[["url_hash"]])
    })


    ##=================================================================================
    ##========================= USER AUTHENTICATION ===================================
    ##=================================================================================

    ## USER <- reactiveValues( logged = FALSE, name="anonymous")
    
    ##================================================================================
    ##====================== INPUT DATA REACTIVE OBJECT ==============================
    ##================================================================================

    currentPGX <- reactiveVal(NULL)
    
    ##inputData <- eventReactive( reload(), {
    inputData <- reactive({
        ## This is wrapper function that loads the ngs object. It also
        ## checks if the user is logged in.
        ##  
        dbg("[LoadingBoard::inputData] ---------- reacted ---------------\n")
        dbg("[LoadingBoard::inputData] AUTHENTICATION=",AUTHENTICATION,"\n")
        dbg("[LoadingBoard::inputData] auth$logged=",auth$logged(),"\n")
        
        ## authenicate user if needed
        ## if(!USER$logged) showLogin()
        ## if(AUTHENTICATION!="none" && USER$logged) showLogin()
        ##if(AUTHENTICATION!="none" && !auth$logged()) return(NULL)
        if(!auth$logged()) return(NULL)

        pgx <- currentPGX()
        dbg("[LoadingBoard::inputData] is.null(pgx)=",is.null(pgx),"\n")
        dbg("[LoadingBoard::inputData] pgx$name = ",pgx$name,"\n")        
        return(pgx)
    })

    observeEvent( input$loadbutton, {

        ## Observe button press
        btn <- isolate(input$loadbutton)
        pgxfile = NULL
        pgxfile = isolate(selectedPGX())
        
        dbg("[LoadingBoard::<loadbutton>] loadbutton=",btn,"\n")
        dbg("[LoadingBoard::<loadbutton>] pgx.selected=",pgxfile)

        if(!is.null(btn) && btn!=0 && !is.null(pgxfile)) {
            ## show loading pop-up
            pgx.showCartoonModal()
        }

        if(is.null(pgxfile) || pgxfile=="" || length(pgxfile)==0) {
            return(NULL)
        }

        message("[LoadingBoard::<loadbutton>] loading pgxfile = ",pgxfile)
        pgx <- loadPGX(pgxfile)

        if(is.null(pgx)) {
            cat("[LoadingBoard::<loadbutton>] ERROR file not found : ",pgxfile,"\n")
            removeModal()
            return(NULL)
        }
        
        ##----------------- update input
        message("[LoadingBoard::<loadbutton>] initializing PGX object")
        pgx <- pgx.initialize(pgx)
        if(is.null(pgx)) {
            cat("[LoadingBoard::<loadbutton>] ERROR in object initialization\n")
            showNotification("ERROR in object initialization!\n")
            removeModal()
            return(NULL)
        }
        if(is.null(pgx$name)) pgx$name <- sub("[.]pgx$","",pgxfile)

        ##----------------- remove modal??
        if(startup_count>0) {
            Sys.sleep(4)
            removeModal()
        }
        
        currentPGX(pgx)
        dbg("[LoadingBoard::<loadbutton>] ready! \n")
    })
    ##}, ignoreNULL=FALSE )
    ##}, ignoreNULL=TRUE )


    ## ================================================================================
    ## =============================== VALUE BOXES UI =================================
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

    getPGXTable <- reactive({ ## reactive 
        
        ## get table of data sets
        ##
        ##
        
        ## Should we read the table from file??
        if(is.null(PGXINFO())) return(NULL)    
        dbg("[LoadingBoard:getPGXTable] *reacted*")
        
        df <- PGXINFO()
        pgxfiles = dir(PGX.DIR, pattern=".pgx$")
        sel <- sub("[.]pgx$","",df$dataset) %in% sub("[.]pgx$","",pgxfiles)
        df <- df[sel,,drop=FALSE]
        dbg("[LoadingBoard:getPGXTable] dim(df)=",dim(df))

        ## Apply filters
        f1=f2=f3=rep(TRUE,nrow(df))
        notnull <- function(x) !is.null(x) && length(x)>0 && x[1]!="" && !is.na(x[1])
        cat("input$flt_datasets = ",input$flt_datasets,"\n")
        cat("input$flt_datatype = ",input$flt_datatype,"\n")
        cat("input$flt_organism = ",input$flt_organism,"\n")
        if(notnull(input$flt_datasets)) f1 <- (df$collection %in% input$flt_datasets)
        if(notnull(input$flt_datatype)) f2 <- (df$datatype %in% input$flt_datatype)
        if(notnull(input$flt_organism)) f3 <- (df$organism %in% input$flt_organism)
        df <- df[which(f1 & f2 & f3),,drop=FALSE]
        
        ##kk = unique(c("dataset","datatype","organism","description",colnames(df)))
        kk = unique(c("dataset","datatype","description","nsamples",
                      "ngenes","nsets","conditions","date",
                      "organism", "collection"))
        kk = intersect(kk,colnames(df))
        df = df[,kk]               
        df = df[order(df$dataset),]   ## sort alphabetically...
        rownames(df) <- NULL
        df
    })


    pgxTable.RENDER <- reactive({

        if(SHOWSPLASH) showStartupModal(once=TRUE)
        
        dbg("[pgxTable.RENDER] reacted")

        df <- getPGXTable()
        req(df)
        dbg("[pgxTable.RENDER] dim(df)=",dim(df))
        
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

    pgxtable_text = "This table contains a general information about all available datasets within the platform. For each dataset, it reports a brief description as well as the total number of samples, genes, gene sets (or pathways), corresponding phenotypes and the collection date."

    pgxtable <- callModule(
        tableModule, id = "pgxtable",
        func = pgxTable.RENDER,
        title = "Datasets",
        height = 640)

    output$pgxtable_UI <- renderUI({    
        fillCol(
            height = 750,
            flex = c(NA,1),
            uiOutput(ns("valueboxes_UI")),
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
    outputOptions(output, "pgxtable_UI", suspendWhenHidden=FALSE) ## important!

    ##================================================================================
    ## Upload data (new)
    ##================================================================================

    output$upload_UI <- renderUI({    
        UploadModuleUI(ns("upload_panel"))
    })

    uploaded_pgx <- UploadModuleServer(
        id = "upload_panel", height = 720,
        ## max.limits = c(samples=20, comparisons=20, genes=8000),
        max.limits = max.limits,
        FILES = FILES        
    )

    observeEvent( uploaded_pgx(), {

        cat("uploaded PGX detected! [LoadingBoard:observe:uploaded_pgx]\n")
        pgx <- uploaded_pgx()
        pgx$collection <- "uploaded"
        ## pgx$owner <- "user"
        currentPGX(pgx)
        selectRows(proxy=dataTableProxy(ns("pgxtable")), selected=NULL)
        
        savedata_button <- NULL
        if(enable_save) {
            
            ##savedata_button <- actionButton(ns("savedata"), "Save my data", icon=icon("save"))
            ##observeEvent( input$savedata, {

            dbg("[LoadingBoard] observeEvent:savedata reacted")        
            ## -------------- save PGX file/object ---------------
            ##pgx <- currentPGX()
            pgx$collection <- "uploaded"
            pgxname <- sub("[.]pgx$","",pgx$name)
            pgxname <- paste0(gsub("[ ]","_",pgxname),".pgx")
            fn  <- file.path(PGX.DIR,pgxname)
            message("[LoadingBoard::@savedata] saving PGX to ",fn)
            message("[LoadingBoard::@savedata] names(pgx)= ",names(pgx),"\n")

            ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ## Note: Currently we use 'ngs' as object name but want to go
            ## towards 'pgx' as standard name. Actually saving as RDS
            ## should be better.
            ngs=pgx
            save(ngs, file=fn)
            remove(ngs)
            
            message("[LoadingBoard::@savedata] updating PGXINFO file")
            new.info <- pgx.updateInfoPGX(PGXINFO(), pgx, remove.old=TRUE)
            Sys.chmod(PGXINFO.FILE, mode="0666")
            write.csv(new.info, file=PGXINFO.FILE)
            PGXINFO(new.info)
            message("[LoadingBoard::@savedata] saved PGXINFO file!")
            
            touchtable(touchtable()+1)
            ##sleep(1)
            ##removeModal()
        }

        
        ## removeModal()
        msg1 <- "<b>Ready!</b><br>You can now start exploring your data. "
        if(enable_save) {
            msg1 <- paste(msg1,"Your data has been saved in your library.")
        }
        showModal( modalDialog(
            HTML(msg1),
            title = NULL,
            size = "s",
            footer = tagList(
                ##savedata_button,
                ## actionButton(ns("sharedata"), "Share with others", icon=icon("share-alt")),
                modalButton("Start!")
            )
        ))

        updateTabsetPanel(session, "tabs",  selected = "Datasets")

    })

    
        

    ##---------------------------------------------------------------
    ##--------------------- modules for UsersMap --------------------
    ##---------------------------------------------------------------
    
    usersmap.RENDER %<a-% reactive({
        
        require(rworldmap)
        require(RColorBrewer)

        df <- ACCESS.LOG$visitors
        
        ##sPDF <- getMap()  
        ##mapCountryData(sPDF, nameColumnToPlot='continent')

        sPDF <- joinCountryData2Map(
            df,
            joinCode = "ISO2",
            nameJoinColumn = "country_code")
        
        par(mai=c(0,0.4,0.2,1),xaxs="i",yaxs="i")
        mapParams <- mapCountryData(
            sPDF, nameColumnToPlot="count",
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

        u <- ACCESS.LOG
        df <- u$visitors
        rownames(df) <-  df$country_name
        tot.users <- sum(df$count)
        freq <- df$count
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
        cat("[LoadingBoard:forum] parenturl =",parenturl,"\n")
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
    ## Board return object
    ##------------------------------------------------
    res <- list(
        inputData = inputData
        ##inputData = currentPGX,
        ##usermode = reactive({ USERMODE() })
    )
    return(res)
}



