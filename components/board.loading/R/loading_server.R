##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

LoadingBoard <- function(id,
                         pgx_dir,
                         pgx,
                         ##authentication="none",
                         auth,
                         limits = c("samples"=1000,"comparisons"=20,
                                    "genes"=20000, "genesets"=10000,
                                    "datasets"=10),
                         enable_upload = TRUE,
                         enable_delete = TRUE,
                         enable_save = TRUE,
                         enable_userdir = TRUE
                         ##force_reload = reactive(0)
                         )
{
  moduleServer(id, function(input, output, session) 
  {
    ns <- session$ns ## NAMESPACE
    dbg("[LoadingBoard] >>> initializing LoadingBoard...")

    loadedDataset <- shiny::reactiveVal(0)  ## counts/trigger dataset upload

    message("[LoadingBoard] in.shinyproxy = ",in.shinyproxy())    
    message("[LoadingBoard] SHINYPROXY_USERNAME = ",Sys.getenv("SHINYPROXY_USERNAME"))
    message("[LoadingBoard] SHINYPROXY_USERGROUPS = ",Sys.getenv("SHINYPROXY_USERGROUPS"))
    message("[LoadingBoard] pgx_dir = ",pgx_dir)
    
    dbg("[LoadingBoard] getwd = ",getwd())

    if(0) {    
      auth <- NULL   ## shared in module
      if(authentication == "password") {
        auth <- shiny::callModule(
          PasswordAuthenticationModule, "auth",
          credentials.file = "CREDENTIALS")
      } else if(authentication == "firebase") {
        auth <- shiny::callModule(FirebaseAuthenticationModule, "auth")
      } else if(authentication == "shinyproxy") {        
        username <- Sys.getenv("SHINYPROXY_USERNAME")
        ##email <- Sys.getenv("SHINYPROXY_EMAIL")        
        auth <- shiny::callModule(NoAuthenticationModule, "auth",
                                  show_modal=TRUE,
                                  username=username, email=username)
      } else if(authentication == "none2") {        
        auth <- shiny::callModule(NoAuthenticationModule, "auth",
                                  show_modal=FALSE)
      } else {
        ##} else if(authentication == "none") {
        auth <- shiny::callModule(NoAuthenticationModule, "auth",
                                  show_modal=TRUE)
      } 
      dbg("[LoadingBoard] names.auth = ",names(auth))
    }
    
    
    ##-----------------------------------------------------------------------------
    ## Description
    ##-----------------------------------------------------------------------------

    shiny::observeEvent( input$module_info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>Load Dataset</strong>"),
            shiny::HTML(module_infotext),
            easyClose = TRUE, size="l" ))
    })

    module_infotext =paste0(
        "This panel shows the available datasets within the platform. The table reports a brief description as well as the total number of samples, genes, gene sets (or pathways), corresponding phenotypes and the creation date.

<br><br><b>Selecting the dataset:</b> Users can select a dataset in the table. The Dataset info shows the information of the dataset of interest and users can load the data by clicking the 'Load dataset' button.

<br><br><br>
<center><iframe width='560' height='315' src='https://www.youtube.com/embed/elwT6ztt3Fo' title='YouTube video player' frameborder='0' allow='accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe><center>

")
    
    ##-----------------------------------------------------------------------------
    ## User interface
    ##-----------------------------------------------------------------------------

    currentSection <- shiny::reactive({
        cdata <- session$clientData
        sub("section-","",cdata[["url_hash"]])
    })

    output$rowselected <- shiny::reactive({
        !is.null(selectedPGX()) && length(selectedPGX())>0
    })
    shiny::outputOptions(output, "rowselected", suspendWhenHidden=FALSE)
    
    observe({
        df <- getPGXINFO()
        datatypes <- sort(setdiff(df$datatype,c(NA,"")))
        organisms <- sort(setdiff(df$organism,c(NA,"")))        
        shiny::updateCheckboxGroupInput(session, "flt_datatype", choices = datatypes)
        shiny::updateCheckboxGroupInput(session, "flt_organism", choices = organisms)
    })
    
    ##-----------------------------------------------------------------------------
    ## READ initial PGX file info
    ##-----------------------------------------------------------------------------
    
    ## reactive value for updating table
    reload_pgxdir <- shiny::reactiveVal(0)
    
    getPGXDIR <- shiny::reactive({
        reload_pgxdir()  ## force reload

        email="../me@company.com"
        email <- auth$email()
        email <- gsub(".*\\/","",email)
        pdir <- pgx_dir  ## from module input

        ##USERDIR=FALSE
        if(enable_userdir) {
            pdir <- paste0(pdir,"/",email)
            if(!is.null(email) && !is.na(email) && email!="") pdir <- paste0(pdir,'/')
            if(!dir.exists(pdir)) {
                dbg("[LoadingBoard:getPGXDIR] userdir does not exists. creating pdir = ",pdir)
                dir.create(pdir)
                dbg("[LoadingBoard:getPGXDIR] copy example pgx")                
                file.copy(file.path(pgx_dir,"example-data.pgx"),pdir)
            }
        }
        pdir
    })
    
    getPGXINFO <- shiny::reactive({
        req(auth)
        if(!auth$logged()) {
            dbg("[LoadingBoard:getPGXINFO] user not logged in!")            
            return(NULL)
        }
        info <- NULL
        pdir <- getPGXDIR()
        info <- pgx.scanInfoFile(pdir, file="datasets-info.csv", verbose=TRUE )
        if(is.null(info)) {
            aa <- rep(NA,9)
            names(aa) = c("dataset","datatype","description","nsamples",
                          "ngenes","nsets","conditions","organism","date")
            info <- data.frame(rbind(aa))[0,]        
        }        
        info
    })

    getFilteredPGXINFO <- shiny::reactive({ 
        
        ## get the filtered table of pgx datasets
        req(auth)
        if(!auth$logged()) {
            dbg("[LoadingBoard:getFilteredPGXINFO] user not logged in! not showing table!")            
            return(NULL)
        }
        df <- getPGXINFO()
        if(is.null(df)) return(NULL)    

        pgxdir <- getPGXDIR()
        pgxfiles = dir(pgxdir, pattern=".pgx$")
        sel <- sub("[.]pgx$","",df$dataset) %in% sub("[.]pgx$","",pgxfiles)
        df  <- df[sel,,drop=FALSE]

        ## Apply filters
        if(nrow(df)>0) {
            f1=f2=f3=rep(TRUE,nrow(df))
            notnull <- function(x) !is.null(x) && length(x)>0 && x[1]!="" && !is.na(x[1])
            if(notnull(input$flt_datatype)) f2 <- (df$datatype %in% input$flt_datatype)
            if(notnull(input$flt_organism)) f3 <- (df$organism %in% input$flt_organism)
            df <- df[which(f1 & f2 & f3),,drop=FALSE]        
            df$date <- as.Date(df$date, format='%Y-%m-%d')
            ## df$date  <- NULL
            df <- df[order(df$dataset),,drop=FALSE]   ## sort alphabetically...
            ##df <- df[order(df$date,decreasing=FALSE),]
            df <- df[order(df$date,decreasing=TRUE),]
            rownames(df) <- nrow(df):1
        }
            
        ##kk = unique(c("dataset","datatype","organism","description",colnames(df)))
        kk = unique(c("dataset","description","datatype","nsamples",
                      "ngenes","nsets","conditions","date","organism"))
        kk = intersect(kk,colnames(df))
        df = df[,kk,drop=FALSE]               
        df
    })
    
    selectedPGX <- shiny::reactive({
        ##sel <- input$pgxtable_rows_selected
        req(pgxtable)
        sel <- pgxtable$rows_selected()
        if(is.null(sel) || length(sel)==0) return(NULL)
        df <- getFilteredPGXINFO()
        if(is.null(df) || nrow(df)==0) return(NULL)        
        pgxfile <- as.character(df$dataset[sel])
        pgxfile <- paste0(sub("[.]pgx$","",pgxfile),".pgx") ## add/replace .pgx
        pgxfile
    })

    ##=============================================================================
    ##========================== OBSERVE/REACT ====================================
    ##=============================================================================
    ##pgxfile="geiger2016-arginine"

    loadPGX <- function(pgxfile) {
        
        req(auth$logged()) 
        if(!auth$logged()) return(NULL)

        pgxfile <- paste0(sub("[.]pgx$","",pgxfile),".pgx") ## add/replace .pgx         
        pgxdir <- getPGXDIR()
        
        pgx.path <- pgxdir[file.exists(file.path(pgxdir,pgxfile))][1]
        pgxfile1 = file.path(pgx.path,pgxfile)
        pgxfile1
        
        pgx <- NULL
        if(file.exists(pgxfile1)) {
            shiny::withProgress(message="Loading data...", value=0.33, {            
              pgx <- local(get(load(pgxfile1,verbose=0)))  ## override any name
            })
        } else {
            message("[LoadingBoard::loadPGX] ***ERROR*** file not found : ",pgxfile)
            return(NULL)
        }
        if(!is.null(pgx)) {
            pgx$name <- pgxfile
            return(pgx)
        } else {
            message("[LoadingBoard::loadPGX] ERROR loading pgx object")
            return(NULL)
        }
    }
    
    output$downloadpgx <- shiny::downloadHandler(
        ##filename = "userdata.pgx",
        filename = function() {
            selectedPGX()
        },
        content = function(file) {
            pgxfile <- selectedPGX()
            cat("[LoadingBoard::loadPGX] pgxfile = ",pgxfile,"\n")
            if(is.null(pgxfile) || pgxfile=="" || length(pgxfile)==0) return(NULL)
            pgx <- loadPGX(pgxfile)
            temp <- tempfile()
            cat("[LoadingBoard::loadPGX] temp = ",temp)            
            save(pgx, file=temp)
            file.copy(temp,file)
        }
    )

    output$downloadzip <- shiny::downloadHandler(
        ##filename = "userdata.zip",
        filename = function() {
            sub("pgx$","zip",selectedPGX())
        },
        content = function(file) {
            pgxfile <- selectedPGX()
            cat("[LoadingBoard::downloadZIP] pgxfile = ",pgxfile,"\n")            
            if(is.null(pgxfile) || pgxfile=="" || length(pgxfile)==0) return(NULL)
            pgxname <- sub("[.]pgx$","",pgxfile)
            pgx <- loadPGX(pgxfile)
            dir.create(tmp <- tempfile())
            tmp2 <- file.path(tmp,pgxname)
            dir.create(tmp2)

            exp.matrix <- sign(pgx$model.parameters$exp.matrix)
            exp.matrix <- contrastAsLabels(exp.matrix) ## new recommended style
            exp.matrix[is.na(exp.matrix)] <- ""
            
            write.csv(round(pgx$counts,digits=2), file=file.path(tmp2, "counts.csv"))
            write.csv(pgx$samples, file=file.path(tmp2, "samples.csv"))
            write.csv(exp.matrix,  file=file.path(tmp2, "contrasts.csv"))
            write.csv(round(pgx$X,digits=4), file=file.path(tmp2, "normalized.csv"))
          
            zipfile <- tempfile(fileext = ".zip")
            zip::zip(zipfile,
                     files=paste0(pgxname,"/",c("counts.csv","samples.csv","contrasts.csv","normalized.csv")),
                     root=tmp)
            ## zip::zip_list(zipfile)
            cat("[LoadingBoard::downloadZIP] zipfile = ",zipfile)                        
            file.copy(zipfile,file)
            remove(pgx); gc();
        }
    )
    
    shiny::observeEvent( input$deletebutton, {
       
        pgxfile <- selectedPGX()
        if(is.null(pgxfile) || pgxfile=="" || length(pgxfile)==0) return(NULL)

        pgx.path <- getPGXDIR()
        pgxfile1 = file.path(pgx.path,pgxfile)
        pgxfile1
        sel <- NULL

        deletePGX <- function() {
            if(input$confirmdelete) {
                cat(">>> deleting",pgxfile,"\n")
                pgxfile2 <- paste0(pgxfile1,"_")  ## mark as deleted
                file.rename(pgxfile1, pgxfile2)                
                reload_pgxdir(reload_pgxdir()+1)          
            } else {
                cat(">>> deletion cancelled\n")
            }
        }

        not.anonymous <- !is.na(auth$name()) && auth$name()!=""         
        allow.delete <- !not.anonymous        
        message("[LoadingBoard::@deletebutton] current user = ",auth$name()," \n")
        message("[LoadingBoard::@deletebutton] allow.delete = ",allow.delete," \n")
        
        allow.delete = TRUE
        if(!allow.delete) {
            message("[LoadingBoard::@deletebutton] WARNING:: ",pgxfile,
                    " not owned by ",auth$name()," \n")
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
    
    
    ##================================================================================
    ##========================== LOAD DATA FROM LIST =================================
    ##================================================================================
    
    load_react <- reactive({
        btn <- input$loadbutton
        query <- parseQueryString(session$clientData$url_search)
        logged <- isolate(auth$logged()) ## avoid reloading when logout/login            
        (!is.null(btn) || !is.null(query[['pgx']])) && logged
    })
    
    shiny::observeEvent( load_react(), {

        if(!load_react()) {
            return(NULL)
        }

        pgxfile = NULL
        
        ## Observe URL query
        query <- parseQueryString(session$clientData$url_search)        
        if(!is.null(query[['pgx']])) {        
            pgxfile <- query[['pgx']]
            pgxfile <- basename(pgxfile)  ## for security
            pgxfile <- paste0(sub("[.]pgx$","",pgxfile),".pgx") ## add/replace .pgx            
        }

        ## Observe button press (over-rides URL query)
        btn <- input$loadbutton
        if(!is.null(btn) && btn!=0) {        
            pgxfile <- selectedPGX()
        }
        
        ## check if file is there
        if(is.na(pgxfile) || is.null(pgxfile) || pgxfile=="" || length(pgxfile)==0) {
            message("[LoadingBoard@load_react] ERROR file not found : ",pgxfile,"\n")
            return(NULL)
        }

        ## During loading show loading pop-up modal
        pgx.showCartoonModal()

        ##---------------------------------------------------------------------
        ##----------------- Loaded PGX object ---------------------------------
        ##---------------------------------------------------------------------
        
        dbg("[LoadingBoard@load_react] loading pgxfile = ",pgxfile)
        loaded_pgx <- loadPGX(pgxfile)
        dbg("[LoadingBoard@load_react] is.null(pgx) = ",is.null(loaded_pgx))
        
        if(is.null(loaded_pgx)) {
            message("[LoadingBoard@load_react] ERROR loading PGX file ",pgxfile,"\n")
            beepr::beep(10)
            shiny::removeModal()            
            return(NULL)
        }
        
        ##----------------- update PGX object ---------------------------------
        dbg("[LoadingBoard@load_react] initializing PGX object")
        loaded_pgx <- pgx.initialize(loaded_pgx)
        dbg("[LoadingBoard@load_react] initialization done!")        


        if(is.null(loaded_pgx)) {
            cat("[LoadingBoard@load_react] ERROR in object initialization\n")
            beepr::beep(10)
            shiny::showNotification("ERROR in object initialization!\n")
            shiny::removeModal()            
            return(NULL)
        }
        loaded_pgx$name <- sub("[.]pgx$","",pgxfile)  ## always use filename
        
        ##----------------- update input --------------------------------------
        loadedDataset(loadedDataset()+1)   ## notify new data uploaded
        
        ## ***NEW*** update PGX from session
        if(1) {
            dbg("[LoadingBoard@load_react] pgx$name = ",loaded_pgx$name)
            dbg("[LoadingBoard@load_react] tracemem(pgx) = ",tracemem(loaded_pgx))

            ## *** EXPERIMENTAL ***. Copying to pgx list to reactiveValues in
            ## session environment.
            dbg("[LoadingBoard@load_react] **** copying current pgx to session.pgx  ****")        
            for(i in 1:length(loaded_pgx)) {
                pgx[[names(loaded_pgx)[i]]] <- loaded_pgx[[i]]
            }
        }

        ##----------------- remove modal on exit?? -------------------------
        ##Sys.sleep(3)
        ##shiny::removeModal()            
        remove(loaded_pgx)
        gc()
        
    })
    ##}, ignoreNULL=FALSE )
    ##}, ignoreNULL=TRUE )


    ## ================================================================================
    ## =============================== VALUE BOXES UI =================================
    ## ================================================================================

    vbox <- function(value, label) {
        div(
            class = "valuebox w-100 border",
            div(
                class = "valuebox-box bg-primary",
                h1(value, class = "valuebox-value"),
                h5(label, class = "valuebox-label")
            )
        )
    }

    output$valuebox1 <- shiny::renderUI({
        pgx <- getFilteredPGXINFO()
        shiny::req(pgx)
        ndatasets = "..."
        ndatasets <- nrow(pgx)
        vbox( ndatasets, "data sets")     
    })

    output$valuebox2 <- shiny::renderUI({
        pgx <- getFilteredPGXINFO()
        shiny::req(pgx)
        ##dbg("valuebox2:: pgx$nsamples=",pgx$nsamples)
        nsamples <- sum(as.integer(pgx$nsamples),na.rm=TRUE)
        vbox( nsamples, "number of samples")     
    })

    output$valuebox3 <- shiny::renderUI({
        pgx <- getFilteredPGXINFO()
        shiny::req(pgx)
        ##dbg("valuebox3:: pgx$nsamples=",pgx$nsamples)
        nvalues <- sum(as.integer(pgx$nsamples) * (as.integer(pgx$ngenes)
            + as.integer(pgx$nsets)),na.rm=TRUE)
        nvalues1 <- format(nvalues, nsmall=, big.mark=" ")  # 1,000.6
        vbox( nvalues1, "data points") 
    })

    output$navheader <- shiny::renderUI({
      fillRow(
        flex=c(NA,1,NA),
        ##h2(input$nav),
        shiny::div(
          id="navheader-current-section",
          HTML("Load dataset &nbsp;"), 
          shiny::actionLink(
            ns("module_info"), "",
            icon=shiny::icon("info-circle"),
            style="color: #ccc;"
            )
        ),        
        shiny::br(),
        shiny::div(selectedPGX(), id="navheader-current-dataset")
      )
    })
    
    ##================================================================================
    ## Data sets
    ##================================================================================
    
    ## reactive value for updating table
    touchtable <- shiny::reactiveVal(0)

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

    
    pgxTable_data <- shiny::reactive({
        
        dbg("[pgxTable.RENDER] reacted")

        ##touchtable()  ## explicit reactive on this
        reload_pgxdir()
        
        df <- getFilteredPGXINFO()
        shiny::req(df)
        dbg("[pgxTable.RENDER] dim(df)=",dim(df))
        
        df$dataset  <- gsub("[.]pgx$"," ",df$dataset)
        df$conditions  <- gsub("[,]"," ",df$conditions)
        df$conditions  <- sapply(as.character(df$conditions), andothers, split=" ", n=5)
        df$description <- shortstring(as.character(df$description),200)
        df$nsets <- NULL
        df$organism <- NULL
        df
    })

    pgxTable_DT <- function() {
        df <- pgxTable_data()
        req(df)
        
        target1 <- grep("date",colnames(df))
        target2 <- grep("description",colnames(df))
        target3 <- grep("conditions",colnames(df))
        target4 <- grep("dataset",colnames(df))        
        
        DT::datatable(
          df,
          class = 'compact cell-border stripe hover',
          rownames=TRUE,
          extensions = c('Scroller'),
          selection = list(mode='single', target='row', selected=1),
          fillContainer = TRUE,
          options=list(
            dom = 'Blfrtip',
            ##columnDefs = list(list(searchable = FALSE, targets = 1)),
            pageLength = 30,
            lengthMenu = c(20, 30, 40, 100),
            scrollX = FALSE,
            ##scrollY =400, ## scroller=TRUE,
            ##scrollY = '100vh', ## scroller=TRUE,
            scrollY = FALSE, 
            deferRender=TRUE,
            autoWidth = TRUE,
            columnDefs = list(
              list(width='60px', targets=target1),
              list(width='30vw', targets=target2)
            )
          )  ## end of options.list 
        )
    }

    pgxTable.RENDER <- function() {
        pgxTable_DT() %>%
            DT::formatStyle(0, target='row', fontSize='13px', lineHeight='95%')        
    }

    pgxTable_modal.RENDER <- function() {
        pgxTable_DT() %>%
            DT::formatStyle(0, target='row', fontSize='20px', lineHeight='95%')        
    }

    info_text = "This table contains a general information about all available datasets within the platform. For each dataset, it reports a brief description as well as the total number of samples, genes, gene sets (or pathways), corresponding phenotypes and the creation date."

    pgxtable <- shiny::callModule(
        tableModule, id = "pgxtable",
        func = pgxTable.RENDER,
        func2 = pgxTable_modal.RENDER,        
        title = "",
        height = c(600,700),
        width = c('100%','100%'),
        info.text = info_text,
        caption2 = info_text
    )
    
    ##------------------------------------------------
    ## Board return object
    ##------------------------------------------------
    res <- list(
        loaded = loadedDataset,
        auth = auth
        ##usermode = shiny::reactive({ USERMODE() })
    )
    return(res)
  })
}
