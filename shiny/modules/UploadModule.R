##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##================================================================================
## Upload data Module
##================================================================================

require(shinyalert)

if(0) {

    OPG = "~/Playground/omicsplayground"
    RDIR = file.path(OPG,"R")
    FILES = file.path(OPG,"lib")
    FILESX = file.path(OPG,"libx")
    PGX.DIR = file.path(OPG,"data")
    source(file.path(RDIR,"pgx-include.R"))  ## pass local vars
    source(file.path(RDIR,"pgx-init.R"))  ## pass local vars

    load(file.path(PGX.DIR,"geiger2016-arginine.pgx"))

    source("UploadModule.R")  ## this file...
    
    pgx <- NULL
    pgx <- gadgetize2(
        UploadModuleUI, UploadModuleServer,
        title = "UploadGadget", height=640, size="l", 
        FILES = "~/Playground/omicsplayground/lib"
    )
    names(pgx)

    pgx <- gadgetize2(
        ComputePgxUI, ComputePgxServer,
        title = "ComputePgxGadget", height=640, size="l", 
        countsRT = reactive(ngs$counts),
        samplesRT = reactive(ngs$samples),
        contrastsRT = reactive(ngs$model.parameters$exp.matrix),
        alertready = TRUE
    )
    pgx <- ComputePgxGadget(ngs$counts, ngs$samples, ngs$model.parameters$exp.matrix)
    names(pgx)

}

UploadModuleUI <- function(id) {
    ns <- NS(id)
    tabsetPanel(
        id = ns("tabs"),
        ## type = "pills",
        ## type = "hidden",
        tabPanel("Upload", uiOutput(ns("upload_UI"))),
        ## tabPanel("Normalize", uiOutput(ns("normalize_UI"))),
        tabPanel("BatchCorrect", uiOutput(ns("batchcorrect_UI"))),
        tabPanel("Contrasts", uiOutput(ns("contrasts_UI"))),
        tabPanel("Compute", uiOutput(ns("compute_UI")))
    )
}

UploadModuleServer <- function(id, height=720, FILES = "../lib", 
                               max.limits = c(samples=100, comparisons=20,
                                              genes=20000, genesets=10000)
                               )
{
    moduleServer(
        id,
        function(input, output, session) {

            message("[UploadModuleServer] called\n")
            ns <- session$ns
            ## ns <- NS(id)
            
            ## Some 'global' reactive variables used in this file
            uploaded <- reactiveValues()

            output$downloadExampleData <- downloadHandler(
                filename = "exampledata.zip",
                content = function(file) {
                    zip = file.path(FILES,"exampledata.zip")
                    file.copy(zip,file)
                }
            )
            
            upload_info = "<h4>User file upload</h4><p>Please prepare the data files in CSV format as listed below. It is important to name the files exactly as shown. The file format must be comma-separated-values (CSV) text. Be sure the dimensions, rownames and column names match for all files. You can download a zip file with example files here: EXAMPLEZIP. You can upload a maximum of <u>LIMITS</u>."
            DLlink = downloadLink(ns("downloadExampleData"),"exampledata.zip")
            upload_info = sub("EXAMPLEZIP", DLlink, upload_info)
            
            ##upload_filetypes = c("text/csv","text/comma-separated-values,text/plain",".csv")
            upload_filetypes = c(".csv",".pgx")                

            limits <- paste(max.limits["samples"],"samples and",
                            max.limits["comparisons"],"comparisons")
            upload_info = sub("LIMITS", limits, upload_info)
            
            ##========================================================================
            ##================================= UI ===================================
            ##========================================================================

            output$upload_UI <- renderUI({
                require(shinyalert)
                useShinyalert()
                fillCol(
                    height = height,
                    ##flex = c(NA,0.05,1.6,0.05,NA),
                    flex = c(NA,0.05,1.6),
                    sidebarLayout(
                        sidebarPanel(
                            width = 3,
                            ## helpText("User file upload"),
                            fileInput(ns("upload_files"), "Choose files",
                                      multiple = TRUE, accept = upload_filetypes),
                            ##checkboxInput(ns("load_example"), "load example data"),
                            shinyWidgets::prettySwitch(ns("load_example"), "load example data"),
                            ##checkboxInput(ns("advanced_mode"),"advanced")
                            ##shinyWidgets::prettySwitch(ns("advanced_mode"),"advanced")
                            shinyWidgets::prettySwitch(ns("advanced_mode"),"batch correction (beta)")
                        ),
                        mainPanel(
                            width = 9,
                            fillRow(
                                flex = c(0.04,1),
                                br(),
                                ##div(HTML(upload_info),style="font-size: 14px;")
                                div(HTML(upload_info))
                            )
                        )
                    ),
                    br(),
                    fillRow(
                        flex = c(0.75,1,1),
                        plotOutput(ns("countStats")) %>% withSpinner(),
                        plotOutput(ns("phenoStats")) %>% withSpinner(),
                        plotOutput(ns("contrastStats")) %>% withSpinner()
                    )
                    ##br(), 
                    ##DT::dataTableOutput(ns("checkTablesOutput"))
                )
            }) ## end-of-renderUI
            
            output$batchcorrect_UI <- renderUI({
                fillCol(
                    height = height, ## width = 1200,
                    BatchCorrectUI(ns("batchcorrect"))
                )                
            })

            output$normalize_UI <- renderUI({
                fillCol(
                    height = height, ## width = 1200,
                    NormalizeCountsUI(ns("normalize"))                    
                )                
            })
            
            output$contrasts_UI <- renderUI({
                fillCol(
                    height = height, ## width = 1200,
                    MakeContrastUI(ns("makecontrast"))
                )                
            })
            outputOptions(output, "contrasts_UI", suspendWhenHidden=FALSE) ## important!!!
            
            output$compute_UI <- renderUI({
                fillCol(
                    height = height, ## width = 1200,
                    ComputePgxUI(ns("compute"))                    
                )                
            })

            ##=====================================================================
            ##============================== TABS =================================
            ##=====================================================================            

            ## !!!!!!!!!!!!! does not work !!!!!!!!!!!!!!!
            observeEvent( uploaded, {
                files.needed = c("counts.csv","samples.csv","contrasts.csv")
                if(all(files.needed %in% names(uploaded))) {
                    showTab("tabs", "Contrasts")
                    showTab("tabs", "Compute")                    
                } else {
                    hideTab("tabs", "Contrasts")
                    hideTab("tabs", "Compute")                    
                }                
            })
            
            ##=====================================================================
            ##========================= OBSERVERS =================================
            ##=====================================================================            
            
            observeEvent( input$advanced_mode, {
                if(input$advanced_mode) {
                    showTab("tabs", "Normalize")   ## NOT YET!!!
                    showTab("tabs", "BatchCorrect")
                } else {
                    hideTab("tabs", "Normalize")
                    hideTab("tabs", "BatchCorrect")
                }
            })

            observeEvent( input$load_example, {
                if(input$load_example) {
                    zipfile = file.path(FILES,"exampledata.zip")
                    readfromzip1 <- function(file) {
                        read.csv(unz(zipfile, file), check.names=FALSE, stringsAsFactors=FALSE,
                                 row.names=1)
                    }
                    readfromzip2 <- function(file) {
                        ## allows for duplicated names
                        df0 <- read.csv(unz(zipfile, file), check.names=FALSE, stringsAsFactors=FALSE)
                        mat <- as.matrix(df0[,-1])
                        rownames(mat) <- as.character(df0[,1])
                        mat
                    }
                    uploaded$counts.csv <- readfromzip2("exampledata/counts.csv")
                    uploaded$samples.csv <- readfromzip1("exampledata/samples.csv")
                    uploaded$contrasts.csv <- readfromzip1("exampledata/contrasts.csv")
                } else {
                    uploaded$counts.csv <- NULL
                    uploaded$samples.csv <- NULL
                    uploaded$contrasts.csv <- NULL
                }
            })
            
            ##=====================================================================
            ##========================= REACTIVES =================================
            ##=====================================================================            

            ##correctedX <- reactive({
            normalized_counts <- NormalizeCountsServerRT(
                id = "normalize",
                counts  = reactive(uploaded$counts.csv),
                height = height
            )
            
            ##correctedX <- reactive({
            correctedX <- BatchCorrectServer(
                id = "batchcorrect",
                X = reactive(uploaded$counts.csv),
                ##X = normalized_counts,  ## NOT YET!!!!
                is.count = TRUE,
                pheno = reactive(uploaded$samples.csv),
                height = height
            )
            
            corrected_counts <- reactive({
                counts <- NULL
                dbg("[UploadModule::corrected_counts] reacted!\n")                
                advanced_mode <- ( length(input$advanced_mode)>0 &&
                                   input$advanced_mode[1]==1 )
                if(advanced_mode) {
                    message("[UploadModule::corrected_counts] using CORRECTED counts\n")
                    out <- correctedX()
                    counts <- pmax(2**out$X-1, 0)
                } else {
                    message("[UploadModule::corrected_counts] using UNCORRECTED counts\n")
                    counts <- uploaded$counts.csv
                }
                counts
            })

            ##mkContrast <- reactive({
            modified_ct <- MakeContrastServerRT(
                id = "makecontrast",
                phenoRT = reactive(uploaded$samples.csv),
                contrRT = reactive(uploaded$contrasts.csv),
                ##countsRT = reactive(uploaded$counts.csv),
                countsRT = corrected_counts,
                height = height
            )
            
            observeEvent( modified_ct(), {
                ## Monitor for changes in the contrast matrix and if
                ## so replace the uploaded reactive values.
                ##
                message("[observe:modified_ct()] reacted...")                
                modct <- modified_ct()
                message("[observe:modified_ct()] dim(modct$contr) = ",dim(modct$contr))
                uploaded$contrasts.csv <- modct$contr
                uploaded$samples.csv   <- modct$pheno

            })

            upload_ok <- reactive({
                check <- checkTables()
                all(check[,"status"]=="OK")
                all(grepl("ERROR",check[,"status"])==FALSE)
            })

            batch_vectors <- reactive({
                correctedX()$B
            })
            
            ##computed_pgx <- ComputePgxServer(
            computed_pgx  <- ComputePgxServer(
                id = "compute",
                ##countsRT = reactive(uploaded$counts.csv),
                countsRT = corrected_counts,
                samplesRT = reactive(uploaded$samples.csv),
                contrastsRT = reactive(uploaded$contrasts.csv),
                batchRT = batch_vectors, 
                enable = upload_ok,
                alertready = FALSE,
                FILES = FILES,
                max.genes = as.integer(max.limits["genes"]),
                max.genesets = as.integer(max.limits["genesets"]),
                height = height
            )

            observeEvent( computed_pgx(), {
                ## Monitor for changes in the computed PGX.
                cat("[UploadModule::computed_pgx] reacted!\n")
                pgx <- computed_pgx()
                if(!is.null(pgx)) {
                    cat("[UploadModule::computed_pgx] updating PGX!\n")
                    uploaded$pgx <- pgx
                }
            })
            
            ##=====================================================================
            ##============================= PLOTS =================================
            ##=====================================================================            

            output$countStats <- renderPlot({

                message("[countStats] renderPlot called")
                ##req(uploaded$counts.csv)                

                check <- checkTables()
                status.ok <- check["counts.csv","status"]                
                message("[countStats] status.ok = ",status.ok)

                if(status.ok!="OK") {
                    frame()
                    status.ds <- check["counts.csv","description"]
                    msg <- paste(toupper(status.ok),"\n\n","(Required) Upload 'counts.csv'",
                                 tolower(status.ds))
                    text(0.5,0.5,paste(strwrap(msg,30),collapse="\n"),col="grey25")
                    box(lty=2, col="grey60")
                    return(NULL)
                }
                
                counts <- uploaded[["counts.csv"]]
                xx <- log2(1 + counts)
                if(nrow(xx)>1000) xx <- xx[sample(1:nrow(xx),1000),,drop=FALSE]
                dc <- melt(xx)
                dc$value[dc$value==0] <- NA
                tt2 <- paste(nrow(counts),"genes x",ncol(counts),"samples")
                ggplot(dc, aes(x=value, color=Var2)) +
                    geom_density() + xlab("log2(1+counts)") +
                    theme( legend.position = "none") +
                    ggtitle("COUNTS", subtitle=tt2)
            })

            output$phenoStats <- renderPlot({

                message("[phenoStats] renderPlot called \n")
                ##req(uploaded$samples.csv)                
                
                check <- checkTables()
                status.ok <- check["samples.csv","status"]                
                if(status.ok!="OK") {
                    frame()
                    status.ds <- check["samples.csv","description"]
                    msg <- paste(toupper(status.ok),"\n\n","(Required) Upload 'samples.csv'",
                                 tolower(status.ds))
                    text(0.5,0.5,paste(strwrap(msg,30),collapse="\n"),col="grey25")
                    box(lty=2, col="grey60")
                    return(NULL)
                }
                
                pheno <- uploaded[["samples.csv"]]                
                px <- head(colnames(pheno),20)  ## show maximum??
                require(inspectdf)
                df <- type.convert(pheno[,px,drop=FALSE])
                vt <- df %>% inspect_types()
                vt
                
                ## discretized continuous variable into 10 bins
                ii <- unlist(vt$col_name[c("numeric","integer")])
                ii
                if(!is.null(ii) && length(ii)) {
                    cat("[UploadModule::phenoStats] discretizing variables:",ii,"\n")
                    df[,ii] <- apply(df[,ii,drop=FALSE], 2, function(x) {
                        if(any(is.infinite(x))) x[which(is.infinite(x))] <- NA
                        cut(x, breaks=10)
                    })
                }
                
                p1 <- df %>% inspect_cat() %>% show_plot()
                tt2 <- paste(nrow(pheno),"samples x",ncol(pheno),"phenotypes")
                ## tt2 <- paste(ncol(pheno),"phenotypes")
                p1 <- p1 + ggtitle("PHENOTYPES", subtitle=tt2) +
                    theme(
                        ##axis.text.x = element_text(size=8, vjust=+5),
                        axis.text.y = element_text(
                            size = 12,
                            margin = ggplot2::margin(0,0,0,25),
                            hjust = 1)
                    )

                message("[UploadModule::phenoStats] done!")
                
                p1
            })
            
            output$contrastStats <- renderPlot({
                
                ##req(uploaded$contrasts.csv)
                ct <- uploaded$contrasts.csv
                has.contrasts <- !is.null(ct) && NCOL(ct)>0
                check <- checkTables()
                status.ok <- check["contrasts.csv","status"]

                message("[output$contrastStats] status.ok = ",status.ok)
                message("[output$contrastStats] has.contrasts = ",has.contrasts)
                message("[output$contrastStats] dim(uploaded$contrasts.csv) = ",
                        dim(uploaded$contrasts.csv))
                
                if( status.ok!="OK" || !has.contrasts) {
                    frame()
                    status.ds <- check["contrasts.csv","description"]
                    msg <- paste(toupper(status.ok),"\n\n","(Optional) Upload 'contrasts.csv'",
                                 tolower(status.ds))
                    ##text(0.5,0.5,"Please upload contrast file 'contrast.csv' with conditions on rows, contrasts as columns")
                    text(0.5,0.5,paste(strwrap(msg,30),collapse="\n"),col="grey25")
                    box(lty=2, col="grey60")
                    return(NULL)
                }

                message("[output$contrastStats] 2 : ")
                
                contrasts <- uploaded$contrasts.csv

                message("[output$contrastStats] 3 : ")                

                ##contrasts <- sign(contrasts)
                ##df <- contrastAsLabels(contrasts)
                df <- contrasts
                px <- head(colnames(df),20)  ## maximum to show??
                df <- data.frame(df[,px,drop=FALSE],check.names=FALSE)
                tt2 <- paste(nrow(contrasts),"samples x",ncol(contrasts),"contrasts")
                ##tt2 <- paste(ncol(contrasts),"contrasts")

                message("[output$contrastStats] 4 : ")
                
                require(inspectdf)
                p1 <- df %>% inspect_cat() %>% show_plot()                    
                p1 <- p1 + ggtitle("CONTRASTS", subtitle=tt2) +
                    theme(
                        ##axis.text.x = element_text(size=8, vjust=+5),
                        axis.text.y = element_text(size = 12,
                                                   margin = ggplot2::margin(0,0,0,25),
                                                   hjust = 1)
                    )

                message("[output$contrastStats] 5 : ")                
                
                p1
            })
            
            
            ##=====================================================================
            ##========================== OBSERVERS ================================
            ##=====================================================================            
            if(0) {
                fn1='~/Downloads/counts.csv'
                fn2='~/Downloads/samples.csv'
                fn3='~/Downloads/contrasts.csv'
                uploadnames=inputnames=c(fn1,fn2)
                uploadnames=inputnames=c(fn1,fn2,fn3)
                
                counts=fread.csv('~/Projects/goutham-csverror/counts.csv')
                samples=read.csv('~/Projects/goutham-csverror/samples.csv',row.names=1)
                contrasts=read.csv('~/Projects/goutham-csverror/contrasts.csv',row.names=1)
                contrasts1=contrasts
                samples1=samples
            }            
            
            ##------------------------------------------------------------------
            ## Main observer for uploaded data files
            ## ------------------------------------------------------------------
            observeEvent( input$upload_files, {
                
                ## Reads in the data files from the file names, checks and
                ## puts in the reactive values object 'uploaded'.
                ##
                message("[upload_files] >>> reading uploaded files")
                message("[upload_files] upload_files$name=",input$upload_files$name)
                message("[upload_files] upload_files$datapath=",input$upload_files$datapath)
                
                ##for(i in 1:length(uploaded)) uploaded[[i]] <- NULL
                uploaded[["pgx"]] <- NULL
                uploaded[["last_uploaded"]] <- NULL
                
                ## read uploaded files
                has.pgx <- any(grepl("[.]pgx$",input$upload_files$name))
                matlist <- list()
                if(has.pgx) {

                    message("[upload_files] extract matrices from PGX")
                    
                    ## If the user uploaded a PGX file, we extract the matrix
                    ## dimensions from the given PGX/NGS object.
                    ##
                    i <- grep("[.]pgx$",input$upload_files$name)
                    load(input$upload_files$datapath[i])  ## load NGS/PGX
                    matlist[["counts.csv"]] <- ngs$counts
                    matlist[["samples.csv"]] <- type.convert(ngs$samples)
                    ##matlist[["contrasts.csv"]] <- ngs$model.parameters$contr.matrix
                    matlist[["contrasts.csv"]] <- ngs$model.parameters$exp.matrix
                    
                    ##if(is.null(ngs$name)) ngs$name <- sub(".pgx$","",input$upload_files$name[i])
                    ##uploaded[["pgx"]] <- ngs
                    
                } else {

                    ## If the user uploaded CSV files, we read in the data
                    ## from the files.
                    ##
                    message("[upload_files] getting matrices from CSV")

                    ii <- grep("csv$",input$upload_files$name)
                    ii <- grep("sample|count|contrast|expression",
                               input$upload_files$name, ignore.case=TRUE)
                    if(length(ii)==0) return(NULL)
                    
                    inputnames  <- input$upload_files$name[ii]
                    uploadnames <- input$upload_files$datapath[ii]

                    message("[upload_files] inputnames = ",inputnames)
                    message("[upload_files] uploadnames = ",uploadnames)                        
                    if(length(uploadnames)>0) {
                        i=1
                        for(i in 1:length(uploadnames)) {
                            fn1 <- inputnames[i]
                            fn2 <- uploadnames[i]
                            matname <- NULL
                            df <- NULL
                            if(grepl("count",fn1, ignore.case=TRUE)) {
                                message("[upload_files] counts.csv : fn1 = ",fn1)
                                ## allows duplicated rownames
                                df0 <- read.csv2(fn2, check.names=FALSE, stringsAsFactors=FALSE)
                                message("[upload_files] counts.csv : 1 : dim(df0) = ",
                                        paste(dim(df0),collapse='x'))
                                if(nrow(df0)>1 && NCOL(df0)>1) {
                                    message("[upload_files] counts.csv : 2 : dim(df0) = ",
                                            paste(dim(df0),collapse='x'))
                                    df <- as.matrix(df0[,-1])
                                    rownames(df) <- as.character(df0[,1])
                                    matname <- "counts.csv"
                                }
                            } else if(grepl("expression",fn1,ignore.case=TRUE)) {
                                message("[upload_files] expression.csv : fn1 = ",fn1)
                                ## allows duplicated rownames
                                df0 <- read.csv2(fn2, check.names=FALSE, stringsAsFactors=FALSE)
                                if(nrow(df0)>1 && NCOL(df0)>1) {
                                    df <- as.matrix(df0[,-1])
                                    rownames(df) <- as.character(df0[,1])
                                    ## convert expression to pseudo-counts
                                    message("[UploadModule::upload_files] converting expression to counts...")
                                    df <- 2**df
                                    matname <- "counts.csv"
                                }
                            } else if(grepl("sample",fn1,ignore.case=TRUE)) {
                                message("[upload_files] samples.csv : fn1 = ",fn1)
                                df <- read.csv2(fn2, row.names=1, check.names=FALSE,
                                                stringsAsFactors=FALSE)
                                df <- type.convert(df)
                                if(nrow(df)>1 && NCOL(df)>=1) {
                                    matname <- "samples.csv"
                                }
                            } else if(grepl("contrast",fn1,ignore.case=TRUE)) {
                                message("[upload_files] contrasts.csv : fn1 = ",fn1)
                                df <- read.csv2(fn2, row.names=1, check.names=FALSE,
                                                stringsAsFactors=FALSE)
                                if(nrow(df)>1 && NCOL(df)>=1) {
                                    matname <- "contrasts.csv"
                                }
                            }
                            if(!is.null(matname)) {
                                matlist[[matname]] <- df
                            }
                        }
                    }            
                }
                
                message("[upload_files] names(matlist) = ",names(matlist))
                if("counts.csv" %in% names(matlist)) {
                    ## Convert to gene names (need for biological effects)
                    message("[upload_files] converting probe names to symbols")
                    X0 <- matlist[['counts.csv']]
                    pp <- rownames(X0)
                    rownames(X0) <- probe2symbol(pp)
                    sel <- !(rownames(X0) %in% c(NA,'','NA'))
                    X0 <- X0[sel,]
                    xx <- tapply(1:nrow(X0), rownames(X0), function(i) colSums(X0[i,,drop=FALSE]))
                    X0 <- do.call(rbind, xx)
                    matlist[['counts.csv']] <- X0
                }
                
                
                ## put the matrices in the reactive values 'uploaded'        
                files.needed = c("counts.csv","samples.csv","contrasts.csv")
                matlist = matlist[ which(names(matlist) %in% files.needed) ]                
                if(length(matlist)>0) {
                    for(i in 1:length(matlist)) {
                        colnames(matlist[[i]]) <- gsub("[\n\t ]","_",colnames(matlist[[i]]))
                        rownames(matlist[[i]]) <- gsub("[\n\t ]","_",rownames(matlist[[i]]))
                        if(names(matlist)[i] %in% c("counts.csv","contrasts.csv")) {
                            matlist[[i]] <- as.matrix(matlist[[i]])
                        } else {
                            matlist[[i]] <- type.convert(matlist[[i]])
                        }
                        m1 <- names(matlist)[i]
                        message("[upload_files] updating matrix ",m1)                        
                        uploaded[[m1]] <- matlist[[i]]
                    }
                    uploaded[["last_uploaded"]] <- names(matlist)
                }
                
                message("[upload_files] done!\n")
            })
                                    
            checkTables <- reactive({        
                ##
                ##
                ##
                
                ## check dimensions
                status = rep("please upload",3)
                files.needed = c("counts.csv","samples.csv","contrasts.csv")        
                names(status) = files.needed
                files.nrow = rep(NA,3)
                files.ncol = rep(NA,3)
                
                for(i in 1:3) {
                    fn = files.needed[i]
                    upfile <- uploaded[[fn]]
                    if(fn %in% names(uploaded) && !is.null(upfile)) {
                        status[i] = "OK"
                        files.nrow[i] = nrow(upfile)
                        files.ncol[i] = ncol(upfile)
                    }
                }

                has.pgx <- ("pgx" %in% names(uploaded))
                if(has.pgx) has.pgx <- has.pgx && !is.null(uploaded[["pgx"]])
                if(has.pgx==TRUE) {
                    
                    ## Nothing to check. Always OK.            
                    
                } else if(!has.pgx) {
                    
                    ## check rownames of samples.csv
                    if(status["samples.csv"]=="OK" && status["counts.csv"]=="OK") {
                        
                        samples1 = uploaded[["samples.csv"]]
                        counts1 = uploaded[["counts.csv"]]
                        a1 <- mean(rownames(samples1) %in% colnames(counts1))
                        a2 <- mean(samples1[,1] %in% colnames(counts1))
                        
                        if(a2 > a1 && NCOL(samples1)>1 ) {
                            message("[UploadModuleServer] getting sample names from first column\n")
                            rownames(samples1) <- samples1[,1]
                            uploaded[["samples.csv"]] <- samples1[,-1,drop=FALSE]
                        }                        
                    }
                    
                    ## check files: matching dimensions
                    if(status["counts.csv"]=="OK" && status["samples.csv"]=="OK") {
                        nsamples   <- ncol(uploaded[["counts.csv"]])
                        ok.samples <- intersect(rownames(uploaded$samples.csv),
                                                colnames(uploaded$counts.csv))
                        n.ok <- length(ok.samples)                        
                        if(n.ok > 0 && n.ok < nsamples) {
                            ## status["counts.csv"]  = "WARNING: some samples with missing annotation)"
                        }
                        if(n.ok > 0) {
                            uploaded[["samples.csv"]] <- uploaded$samples.csv[ok.samples,,drop=FALSE]
                            uploaded[["counts.csv"]]  <- uploaded$counts.csv[,ok.samples,drop=FALSE]
                        }
                        if(n.ok == 0) {
                            status["counts.csv"]  = "ERROR: colnames do not match (with samples)"
                            status["samples.csv"] = "ERROR: rownames do not match (with counts)"
                        }
                    }
                    
                    if(status["contrasts.csv"]=="OK" && status["samples.csv"]=="OK") {
                        samples1   <- uploaded[["samples.csv"]]
                        contrasts1 <- uploaded[["contrasts.csv"]]
                        old1 = ("group" %in% colnames(samples1) &&
                                nrow(contrasts1) < nrow(samples1) &&
                                all(rownames(contrasts1) %in% samples1$group)
                        )
                        old2 = all(rownames(contrasts1)==rownames(samples1)) &&
                            all(unique(as.vector(contrasts1)) %in% c(-1,0,1,NA))

                        old.style <- (old1 || old2)
                        if(old.style && old1) {

                            message("[UploadModule] WARNING: converting old1 style contrast to new format")
                            new.contrasts <- samples1[,0]
                            if(NCOL(contrasts1)>0) {
                                new.contrasts <- contrastAsLabels(contrasts1)
                                grp = as.character(samples1$group)
                                new.contrasts <- new.contrasts[grp,,drop=FALSE]
                                rownames(new.contrasts) <- rownames(samples1)
                            }
                            message("[UploadModule] old.ct1 = ",paste(contrasts1[,1],collapse=' '))
                            message("[UploadModule] old.nn = ",paste(rownames(contrasts1),collapse=' '))
                            message("[UploadModule] new.ct1 = ",paste(new.contrasts[,1],collapse=' '))
                            message("[UploadModule] new.nn = ",paste(rownames(new.contrasts),collapse=' '))
                            
                            contrasts1 <- new.contrasts
                        }
                        if(old.style && old2 ) {
                            message("[UploadModule] WARNING: converting old2 style contrast to new format")
                            new.contrasts <- samples1[,0]
                            if(NCOL(contrasts1)>0) {
                                new.contrasts <- contrastAsLabels(contrasts1)
                                rownames(new.contrasts) <- rownames(samples1)
                            }
                            contrasts1 <- new.contrasts
                        }

                        message("[UploadModule] 1 : dim.contrasts1 = ",dim(contrasts1))
                        if(NCOL(contrasts1)>0) {
                            ## always clean up
                            contrasts1 <- apply(contrasts1,2,as.character)
                            message("[UploadModule] 2 : dim.contrasts1 = ",dim(contrasts1))
                            rownames(contrasts1) <- rownames(samples1)
                            message("[UploadModule] 3 : dim.contrasts1 = ",dim(contrasts1))
                            for(i in 1:ncol(contrasts1)) {
                                isz = (contrasts1[,i] %in% c(NA,"NA","NA ",""," ","  ","   "," NA"))
                                if(length(isz)) contrasts1[isz,i] <- NA
                            }
                        } else {
                            ## contrasts1 <- NULL
                        }
                        uploaded[["contrasts.csv"]] <- contrasts1                        
                    }

                    MAXSAMPLES   = 25
                    MAXCONTRASTS = 5
                    MAXSAMPLES   = as.integer(max.limits["samples"])
                    MAXCONTRASTS = as.integer(max.limits["comparisons"])

                    ## check files: maximum contrasts allowed
                    if(status["contrasts.csv"]=="OK") {
                        if( ncol(uploaded[["contrasts.csv"]]) > MAXCONTRASTS ) {
                            status["contrasts.csv"] = paste("ERROR: max",MAXCONTRASTS,"contrasts allowed")
                        }
                    }
                    
                    ## check files: maximum samples allowed
                    if(status["counts.csv"]=="OK" && status["samples.csv"]=="OK") {
                        if( ncol(uploaded[["counts.csv"]]) > MAXSAMPLES ) {
                            status["counts.csv"]  = paste("ERROR: max",MAXSAMPLES," samples allowed")
                        }
                        if( nrow(uploaded[["samples.csv"]]) > MAXSAMPLES ) {
                            status["samples.csv"] = paste("ERROR: max",MAXSAMPLES,"samples allowed")
                        }
                    }
                    
                    ## check files: must have group column defined
                    if(status["samples.csv"]=="OK" && status["contrasts.csv"]=="OK") {
                        samples1   = uploaded[["samples.csv"]]
                        contrasts1 = uploaded[["contrasts.csv"]]
                        if(!all(rownames(contrasts1) %in% rownames(samples1))) {
                            status["contrasts.csv"] = "ERROR: contrasts do not match samples"
                        }
                    }                    
                    
                } ## end-if-from-pgx
                
                e1 <- grepl("ERROR",status["samples.csv"])
                e2 <- grepl("ERROR",status["contrasts.csv"]) 
                e3 <- grepl("ERROR",status["counts.csv"])
                s1 <- "samples.csv" %in% uploaded$last_uploaded
                s2 <- "contrasts.csv" %in% uploaded$last_uploaded
                s3 <- "counts.csv" %in% uploaded$last_uploaded
                
                if( e1 || e2 || e3 ) {
                    message("[checkTables] ERROR in samples table : e1 = ",e1)
                    message("[checkTables] ERROR in contrasts table : e2 = ",e2)
                    message("[checkTables] ERROR in counts table : e2 = ",e3)

                    if(e1 && !s1) {
                        uploaded[["samples.csv"]] <- NULL
                        status["samples.csv"] = "please upload"
                    }
                    if(e2 && !s2) {
                        uploaded[["contrasts.csv"]] <- NULL
                        status["contrasts.csv"] = "please upload"
                    }
                    if(e3 && !s3) {
                        uploaded[["counts.csv"]] <- NULL
                        status["counts.csv"] = "please upload"
                    }
                }


                if( !is.null(uploaded$contrasts.csv) &&
                    (is.null(uploaded$counts.csv) ||
                     is.null(uploaded$samples.csv)) )
                {
                    uploaded[["contrasts.csv"]] <- NULL
                    status["contrasts.csv"] = "please upload"
                }
                
                
                ## check files
                description = c(
                    "Count/expression file with gene on rows, samples as columns",
                    "Samples file with samples on rows, phenotypes as columns",
                    ## "Gene information file with genes on rows, gene info as columns.",
                    "Contrast file with conditions on rows, contrasts as columns"        
                )
                df <- data.frame(
                    filename = files.needed,
                    description = description,
                    nrow = files.nrow,
                    ncol = files.ncol,
                    status = status
                )
                rownames(df) <- files.needed
                
                ## deselect
                ## selectRows(proxy = dataTableProxy("pgxtable"), selected=NULL)
                return(df)    
            })
            
            output$checkTablesOutput <- DT::renderDataTable({
                ## Render the upload status table
                ##
                if(!input$advanced_mode) return(NULL)
                df <- checkTables()
                dt <- datatable(
                    df,
                    rownames=FALSE,
                    selection = 'none',
                    class="compact cell-border",
                    options = list(
                        dom = 't'
                    )                                
                ) %>%
                    DT::formatStyle(0, target='row', fontSize='12px', lineHeight='100%')                
            })               

            ##========================================================================
            ## return results as reactive object
            ##========================================================================            
            ## return(reactive(uploaded$pgx))  ## pointing to reactive results object
            return(computed_pgx)
            
        }  ## end function
    )  ## end moduleServer
} ## end UploadModuleServer
