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
    load(file.path(PGX.DIR,"carr2019-jev.pgx"))
    load(file.path(PGX.DIR,"GSE10846-dlbcl-nc.pgx"))
    load(file.path(PGX.DIR,"geiger2016-arginine.pgx"))

    source(file.path(RDIR,"pgx-init.R"))  ## pass local vars
    source(file.path(RDIR,"pgx-include.R"))  ## pass local vars
    source(file.path(RDIR,"pgx-ui.R"))  ## pass local vars  
    ##source("../../R/pgx-init.R")      ## pass local vars
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
        tabPanel("Normalize", uiOutput(ns("normalize_UI"))),
        tabPanel("BatchCorrect", uiOutput(ns("batchcorrect_UI"))),
        tabPanel("Contrasts", uiOutput(ns("contrasts_UI"))),
        tabPanel("Compute", uiOutput(ns("compute_UI")))
    )
}

UploadModuleServer <- function(id, height=720, FILES = "../lib", 
                               max.limits = c(samples=20, comparisons=20, genes=8000)
                               )
{
    moduleServer(
        id,
        function(input, output, session) {

            cat("[UploadModuleServer] called\n")
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
                            checkboxInput(ns("load_example"), "load example data"),
                            ##checkboxInput(ns("advanced_mode"),"advanced")
                            shinyWidgets::prettySwitch(ns("advanced_mode"),"advanced")
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
                        plotOutput(ns("countStats"))  %>% withSpinner(),
                        plotOutput(ns("phenoStats"))  %>% withSpinner(),
                        plotOutput(ns("contrastStats"))  %>% withSpinner()
                    )
                    ##br(), 
                    ##DT::dataTableOutput(ns("statusTableOutput"))
                )
            }) ## end-of-renderUI
            
            output$batchcorrect_UI <- renderUI({
                fillCol(
                    height = height, ## width = 1200,
                    SuperBatchCorrectUI(ns("batchcorrect"))
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
            
            output$compute_UI <- renderUI({
                fillCol(
                    height = height, ## width = 1200,
                    ComputePgxUI(ns("compute"))                    
                )                
            })

            ##=====================================================================
            ##============================== TABS =================================
            ##=====================================================================            

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
                    ## showTab("tabs", "Normalize")   ## NOT YET!!!
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
            correctedX <- SuperBatchCorrectServer(
                id = "batchcorrect",
                X = reactive(uploaded$counts.csv),
                ## X = normalized_counts,  ## NOT YET!!!!
                is.count = TRUE,
                ##pheno = reactive(uploaded$samples.csv),
                pheno = reactive(uploaded$samples.csv),
                height = height
            )
            
            ##mkContrast <- reactive({
            modified_ct <- MakeContrastServerRT(
                id = "makecontrast",
                phenoRT = reactive(uploaded$samples.csv),
                contrRT = reactive(uploaded$contrasts.csv),
                height = height
            )

            observeEvent( modified_ct(), {
                ## Monitor for changes in the contrast matrix and if
                ## so replace the uploaded reactive values.
                ##
                cat("[UploadModule::modified_ct] reacted!\n")
                ct <- modified_ct()
                if(!is.null(ct)) {
                    cat("[UploadModule::modified_ct] updating contrast!\n")
                    cat("[UploadModule::modified_ct] dim(ct)=",dim(ct),"\n")
                    uploaded$contrasts.csv <- ct
                }
            })

            upload_ok <- reactive({
                statusTab <- statusTable()
                all(statusTab[,"status"]=="OK")
            })

            corrected_counts <- reactive({
                counts <- NULL
                if(input$advanced_mode) {
                    cat("[UploadModule::corrected_counts] using CORRECTED counts\n")
                    out <- correctedX()
                    counts <- pmax(2**out$X -1, 0)
                } else {
                    cat("[UploadModule::corrected_counts] using UNCORRECTED counts\n")
                    counts <- uploaded$counts.csv
                }
                counts
            })
                        
            ##computed_pgx <- ComputePgxServer(
            computed_pgx  <- ComputePgxServer(
                id = "compute",
                ##countsRT = reactive(uploaded$counts.csv),
                countsRT = corrected_counts,
                samplesRT = reactive(uploaded$samples.csv),
                contrastsRT = reactive(uploaded$contrasts.csv),
                enable = upload_ok,
                alertready = FALSE,
                FILES = FILES,
                max.features = as.integer(max.limits["genes"]),
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

                cat("[countStats] renderPlot called \n")
                ##req(uploaded$counts.csv)                

                statusTab <- statusTable()
                status.ok <- statusTab["counts.csv","status"]                
                if(status.ok!="OK") {
                    frame()
                    status.ds <- statusTab["counts.csv","description"]
                    msg <- paste(toupper(status.ok),"\n\n","Please upload 'counts.csv'",
                                 tolower(status.ds))
                    ##text(0.5,0.5,"Please upload contrast file 'contrast.csv' with conditions on rows, contrasts as columns")
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

                cat("[phenoStats] renderPlot called \n")
                ##req(uploaded$samples.csv)                
                
                statusTab <- statusTable()
                status.ok <- statusTab["samples.csv","status"]                
                if(status.ok!="OK") {
                    frame()
                    status.ds <- statusTab["samples.csv","description"]
                    msg <- paste(toupper(status.ok),"\n\n","Please upload 'samples.csv'",
                                 tolower(status.ds))
                    ##text(0.5,0.5,"Please upload contrast file 'contrast.csv' with conditions on rows, contrasts as columns")
                    text(0.5,0.5,paste(strwrap(msg,30),collapse="\n"),col="grey25")
                    box(lty=2, col="grey60")
                    return(NULL)
                }
                
                pheno <- uploaded[["samples.csv"]]
                px <- head(colnames(pheno),20)  ## maximum??
                require(inspectdf)
                df <- type.convert(pheno[,px,drop=FALSE])
                vt <- df %>% inspect_types()
                vt

                ## discretized continuous variable into 10 bins
                ii <- unlist(vt$col_name[c("numeric","integer")])
                ii
                if(!is.null(ii) && length(ii)) {
                    cat("[UploadModule::phenoStats] discretizing variables:",colnames(df)[ii],"\n")
                    df[,ii] <- apply(df[,ii,drop=FALSE], 2, function(x) cut(x, breaks=10))
                }

                cat("[UploadModule::phenoStats] dim(df)=",dim(df),"\n")
                cat("[UploadModule::phenoStats] class(df)=",class(df),"\n")
                
                p1 <- df %>% inspect_cat() %>% show_plot()

                cat("[UploadModule::phenoStats] class(p1)=",class(p1),"\n")

                tt2 <- paste(nrow(pheno),"samples x",ncol(pheno),"phenotypes")
                ## tt2 <- paste(ncol(pheno),"phenotypes")
                p1 <- p1 + ggtitle("PHENOTYPES", subtitle=tt2) +
                    theme(
                        ##axis.text.x = element_text(size=8, vjust=+5),
                        axis.text.y = element_text(size = 12,
                                                   margin = margin(0,0,0,25),
                                                   hjust = 1)
                    )
                p1
            })
            
            output$contrastStats <- renderPlot({
                
                has.contrasts <- !is.null(uploaded$contrasts.csv) &&
                    NCOL(uploaded$contrasts.csv)>0

                statusTab <- statusTable()
                status.ok <- statusTab["contrasts.csv","status"]
                
                if(status.ok!="OK") {
                    frame()
                    status.ds <- statusTab["contrasts.csv","description"]
                    msg <- paste(toupper(status.ok),"\n\n","Please upload 'contrasts.csv'",
                                 tolower(status.ds))
                    ##text(0.5,0.5,"Please upload contrast file 'contrast.csv' with conditions on rows, contrasts as columns")
                    text(0.5,0.5,paste(strwrap(msg,30),collapse="\n"),col="grey25")
                    box(lty=2, col="grey60")
                    return(NULL)
                }
                
                require(inspectdf)
                req(uploaded$contrasts.csv, uploaded$samples.csv)
                pheno <- uploaded$samples.csv
                contrasts <- uploaded$contrasts.csv
                
                ## convert to experiment matrix if not already...
                if(!all(rownames(contrasts)==rownames(pheno))) {
                    contrasts <- pgx.expMatrix(pheno, contrasts)
                }
                ##contrasts <- sign(contrasts)
                df <- contrastAsLabels(contrasts)
                px <- head(colnames(df),20)  ## maximum to show??
                df <- data.frame(df[,px,drop=FALSE],check.names=FALSE)
                tt2 <- paste(nrow(contrasts),"samples x",ncol(contrasts),"contrasts")
                ##tt2 <- paste(ncol(contrasts),"contrasts")

                cat("[UploadModule::phenoStats] dim(df)=",dim(df),"\n")
                cat("[UploadModule::phenoStats] tt2=",tt2,"\n")

                p1 <- df %>% inspect_cat() %>% show_plot()                    

                cat("[UploadModule::phenoStats] class(p1)=",class(p1),"\n")
                
                p1 <- p1 + ggtitle("CONTRASTS", subtitle=tt2) +
                    theme(
                        ##axis.text.x = element_text(size=8, vjust=+5),
                        axis.text.y = element_text(size = 12,
                                                   margin = margin(0,0,0,25),
                                                   hjust = 1)
                    )
                p1

            })
            
            
            ##=====================================================================
            ##========================== OBSERVERS ================================
            ##=====================================================================            

            observeEvent( input$autocontrast, {
                req(uploaded$samples.csv)

                cat("[!autocontrast] >>> creating autocontrast\n")                
                df <- uploaded$samples.csv
                cat("[!autocontrast] dim(df)=",dim(df),"\n")                
                ct <- pgx.makeAutoContrast(
                    df, mingrp=3, slen=20, ref=NULL, fix.degenerate=FALSE)
                rownames(ct$exp.matrix) <- rownames(df)

                cat("[!autocontrast] updating contrasts...\n")                                
                uploaded[["contrasts.csv"]] <- ct$exp.matrix
            })
            
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
                
                ## read uploaded files
                has.pgx <- any(grepl("[.]pgx$",input$upload_files$name))
                matlist <- list()
                if(has.pgx) {
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
                    ii <- grep("csv$",input$upload_files$name)
                    inputnames  <- input$upload_files$name[ii]
                    uploadnames <- input$upload_files$datapath[ii]
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
                                df <- read.csv(fn2, row.names=1, check.names=FALSE,
                                               stringsAsFactors=FALSE)
                            }
                            matlist[[ inputnames[i] ]] <- df
                        }
                    }            
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
                        uploaded[[names(matlist)[i]]] <- matlist[[i]]
                    }
                }

                cat("[upload_files] dim(counts.csv)=",dim(uploaded$counts.csv),"\n")
                cat("[upload_files] dim(samples.csv)=",dim(uploaded$samples.csv),"\n")
                cat("[upload_files] dim(contrasts.csv)=",dim(uploaded$contrasts.csv),"\n")
                
                cat("[upload_files] done!\n")
            })
                                    
            statusTable <- reactive({        
                
                message("[statusTable] ********** CALLED ***************")

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
                    
                    ## check files: matching dimensions
                    if(status["counts.csv"]=="OK" && status["samples.csv"]=="OK") {
                        if(!all( sort(colnames(uploaded[["counts.csv"]])) ==
                                 sort(rownames(uploaded[["samples.csv"]])) )) {
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
                        has.group <- "group" %in% colnames(samples1)
                        matching.group <- all(rownames(contrasts1) %in% samples1$group)
                        matching.samples <- all(rownames(contrasts1) %in% rownames(samples1))
                        ##if(!has.group) {
                        ##    status["samples.csv"] = "ERROR: missing 'group' column"
                        ##}
                        if(!matching.group && !matching.samples) {
                            status["contrasts.csv"] = "ERROR: contrasts do not match samples"
                        }
                    }
                    
                } ## end-if-from-pgx
                
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
            
            output$statusTableOutput <- DT::renderDataTable({
                ## Render the upload status table
                ##
                if(!input$advanced_mode) return(NULL)
                df <- statusTable()
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
