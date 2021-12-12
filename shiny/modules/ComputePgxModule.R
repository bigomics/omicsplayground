##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

if(0) {

    source("~/Playground/omicsplayground/R/pgx-include.R")
    load("~/Playground/omicsplayground/data/GSE10846-dlbcl-nc.pgx")    

    PgxComputeGadget(X=ngs$X, pheno=ngs$samples)

    out <- gadgetize2(
        ComputePgxUI, ComputePgxServer,
        title = "UploadGadget", height=640, size="l", 
        X = ngs$X, pheno=ngs$samples )
    names(out)
    
}

ComputePgxGadget <- function(counts, samples, contrasts, height=720) {
    gadgetize(
        ComputePgxUI, ComputePgxServer,
        title = "ComputePGX",
        countsRT = shiny::reactive(counts),
        samplesRT = shiny::reactive(samples),
        contrastsRT = shiny::reactive(contrasts)
    )
}

ComputePgxUI <- function(id) {
    ns <- shiny::NS(id)
    shiny::uiOutput(ns("UI"))
}

ComputePgxServer <- function(id, countsRT, samplesRT, contrastsRT, batchRT, metaRT,
                             FILES, pgx.dirRT, enable_button = TRUE, alertready = TRUE, 
                             max.genes = 20000, max.genesets = 10000,
                             max.datasets = 100, height = 720 )
{
    shiny::moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns

            dbg("ComputePgxServer init")

            ## statistical method for GENE level testing
            GENETEST.METHODS = c("ttest","ttest.welch","voom.limma","trend.limma","notrend.limma",
                                 "deseq2.wald","deseq2.lrt","edger.qlf","edger.lrt")
            GENETEST.SELECTED = c("trend.limma","deseq2.wald","edger.qlf")
            
            ## statistical method for GENESET level testing
            GENESET.METHODS = c("fisher","ssgsea","gsva", "spearman", "camera", "fry",
                                ##"plage","enricher","gsea.permPH","gsea.permGS","gseaPR",
                                "fgsea")
            GENESET.SELECTED = c("fisher","gsva","fgsea")

            ## batch correction and extrs methods
            EXTRA.METHODS = c("deconv", "drugs", "wordcloud","connectivity")
            EXTRA.NAMES = c("celltype deconvolution", "drugs connectivity",
                            "wordcloud","experiment similarity")
            EXTRA.SELECTED = c("deconv","drugs","wordcloud","connectivity")

            DEV.METHODS = c("noLM.prune")
            DEV.NAMES = c("noLM + prune")
            DEV.SELECTED = c()

            output$UI <- shiny::renderUI({
                shiny::fillCol(
                    height = height,
                    flex = c(0.2,NA,0.05,1.3),
                    shiny::br(),
                    shiny::fluidRow(
                        shiny::column(
                            12, align="center", offset=0,
                            shiny::tags$table(
                                     style="width:100%;vertical-align:top;padding:4px;",
                                     shiny::tags$tr(
                                              shiny::tags$td("", width="300"),
                                              shiny::tags$td("Name", width="100"),
                                              shiny::tags$td(shiny::textInput(
                                                       ns("upload_name"),NULL, ##"Dataset:",
                                                       placeholder="Name of your dataset"),
                                                      width="600"
                                                      ),
                                              shiny::tags$td("", width="120")
                                          ),
                                     shiny::tags$tr(
                                              shiny::tags$td(""),
                                              shiny::tags$td("Datatype"),
                                              shiny::tags$td(shiny::selectInput(
                                                       ns("upload_datatype"), NULL,
                                                       choices = c("RNA-seq","scRNA-seq", "proteomics",
                                                                   "mRNA microarray","other"))
                                                      ),
                                              shiny::tags$td("")
                                          ),
                                     shiny::tags$tr(
                                              shiny::tags$td(""),
                                              shiny::tags$td("Description"),
                                              shiny::tags$td(shiny::div(shiny::textAreaInput(
                                                       ns("upload_description"), NULL, ## "Description:",
                                                       placeholder="Give a short description of your dataset",
                                                       height=100, resize='none'),
                                                       style="margin-left: 0px;"
                                                       )),
                                              shiny::tags$td("")
                                          )
                                 ),
                            shiny::br(),
                            shiny::div(
                                shiny::actionButton(ns("compute"),"Compute!",icon=icon("running"),
                                                    class="run-button"),
                                shiny::br(),br(),
                                shiny::actionLink(ns("options"), "Computation options",
                                                  icon=icon("cog", lib="glyphicon")),
                                style = "padding-right: 80px;"
                                )
                        )
                    ),
                    shiny::br(),
                    shiny::conditionalPanel(
                        "input.options%2 == 1", ns=ns,
                        shiny::fillRow(
                            shiny::wellPanel(
                                shiny::checkboxGroupInput(
                                    ns('filter_methods'),
                                    'Feature filtering:',
                                    choiceValues =
                                        c("only.hugo",
                                          "only.proteincoding",
                                          "remove.unknown",
                                          "remove.notexpressed"
                                          ## "excl.immuno"
                                          ## "excl.xy"
                                          ),
                                    choiceNames =
                                        c("convert to HUGO",
                                          "protein-coding only",
                                          "remove Rik/ORF/LOC genes",
                                          "remove not-expressed"
                                          ##"Exclude immunogenes",
                                          ##"Exclude X/Y genes"
                                          ),
                                    selected = c(
                                        "only.hugo",
                                        "only.proteincoding",
                                        "remove.unknown",
                                        "remove.notexpressed"
                                    )
                                )
                            ),
                            shiny::wellPanel(
                                shiny::checkboxGroupInput(
                                    ns('gene_methods'),
                                    'Gene tests:',
                                    GENETEST.METHODS,
                                    selected = GENETEST.SELECTED
                                )
                            ),
                            shiny::wellPanel(
                                shiny::checkboxGroupInput(
                                    ns('gset_methods'),
                                    'Enrichment methods:',
                                    GENESET.METHODS,
                                    selected = GENESET.SELECTED
                                ),
                            ),
                            shiny::wellPanel(
                                shiny::checkboxGroupInput(
                                    ns('extra_methods'),
                                    'Extra analysis:',
                                    choiceValues = EXTRA.METHODS,
                                    choiceNames = EXTRA.NAMES,
                                    selected = EXTRA.SELECTED
                                )
                            ),
                            shiny::wellPanel(
                                shiny::checkboxGroupInput(
                                    ns('dev_options'),
                                    'Developer options:',
                                    choiceValues = DEV.METHODS,
                                    choiceNames = DEV.NAMES,
                                    selected = DEV.SELECTED
                                )
                            )
                        ) ## end of fillRow
                    ) ## end of conditional panel
                ) ## end of fill Col
            })
            shiny::outputOptions(output, "UI", suspendWhenHidden=FALSE) ## important!!!            
            

            if(FALSE) {
                shiny::observeEvent( input$gene_methods, {
                    if(length(input$gene_methods) > 3){
                        shiny::updateCheckboxGroupInput(session, "gene_methods",
                                                        selected= tail(input$gene_methods,3))
                    }
                    if(length(input$gene_methods) < 1){
                        shiny::updateCheckboxGroupInput(session, "gene_methods", selected= "ttest")
                    }
                })
            }
            
            shiny::observeEvent( input$options, {
                ## shinyjs::disable(ns("gene_methods2"))
            })
            
            shiny::observeEvent( enable_button(), {
                ## NEED CHECK. not working... 
                ##
                if(!enable_button()){
                    message("[ComputePgxServer:@enable] disabling compute button")
                    shinyjs::disable(ns("compute"))
                } else {
                    message("[ComputePgxServer:@enable] enabling compute button")
                    shinyjs::enable(ns("compute"))
                }
            })
            
            shiny::observeEvent( metaRT(), {
                
                dbg("[ComputePgxServer:@metaRT] parsing meta information...")
                meta <- metaRT()

                dbg("[ComputePgxServer:@metaRT] names.meta =",names(meta))
                
                if(!is.null(meta[['name']])) {
                    dbg("[ComputePgxServer:@metaRT] NS.upload_name = ",ns("upload_name"))
                    dbg("[ComputePgxServer:@metaRT] meta.name => ",meta[['name']])
                    shiny::updateTextInput(session, "upload_name", value=meta[['name']])
                    ## shiny::updateTextInput(session, ns("upload_name"), value=meta[['name']])                    
                }
                if(!is.null(meta[['description']])) {
                    dbg("[ComputePgxServer:@metaRT] NS.upload_description = ",ns("upload_description"))
                    dbg("[ComputePgxServer:@metaRT] meta.description => '",meta[['description']],"'")
                    shiny::updateTextAreaInput(session, "upload_description", value=meta[['description']])
                    ##shiny::updateTextAreaInput(session, ns("upload_description"), value=meta[['description']])                    
                }
                
            })

            ##------------------------------------------------------------------
            ## After confirmation is received, start computing the PGX
            ## object from the uploaded files
            ## ------------------------------------------------------------------
            computedPGX  <- shiny::reactiveVal(NULL)
            
            shiny::observeEvent( input$compute, {
                
                message("[ComputePgxServer::@compute] reacted!")
                ## shiny::req(input$upload_hugo,input$upload_filtergenes)
                
                ##-----------------------------------------------------------
                ## Check validity 
                ##-----------------------------------------------------------                
                if(!enable_button()) {
                    message("[ComputePgxServer:@compute] WARNING:: *** NOT ENABLED ***")
                    return(NULL)
                }

                pgxdir <- pgx.dirRT()
                dbg("[ComputePgxServer::@compute] pgxdir  = ", pgxdir )

                numpgx <- length(dir(pgxdir, pattern="*.pgx$"))
                dbg("[ComputePgxServer::@compute] numpgx  = ", numpgx )                
               
                if(numpgx >= max.datasets) {
                    msg = "Your storage is full. You have NUMPGX pgx files in your data folder and your quota is LIMIT datasets. Please delete some datasets or consider buying extra storage."
                    msg <- sub("NUMPGX",numpgx,msg)
                    msg <- sub("LIMIT",max.datasets,msg)
                    shinyalert::shinyalert("WARNING",msg)
                    return(NULL)
                }

                has.contrasts <- !is.null(contrastsRT()) && NCOL(as.matrix(contrastsRT()))>0
                if(!has.contrasts) {
                    shinyalert::shinyalert("ERROR","You must define at least 1 contrast")
                    return(NULL)
                }

                has.name <- input$upload_name != ""
                has.description <- input$upload_description != ""
                if(!has.name || !has.description) {
                    shinyalert::shinyalert("ERROR","You must give a dataset name and description")
                    return(NULL)
                }
                
                ##-----------------------------------------------------------
                ## Retrieve the most recent matrices from reactive values
                ##-----------------------------------------------------------                
                counts    <- countsRT()
                samples   <- samplesRT()
                samples   <- data.frame(samples, stringsAsFactors=FALSE, check.names=FALSE)
                contrasts <- as.matrix(contrastsRT())
                
                dbg("[ComputePgxServer:@enable] ct1 = ", contrasts[,1])
                
                ## contrasts[is.na(contrasts)] <- 0
                ## contrasts[is.na(contrasts)] <- ""                
                ##!!!!!!!!!!!!!! This is blocking the computation !!!!!!!!!!!
                ##batch  <- batchRT() ## batch correction vectors for GLM
                
                ##-----------------------------------------------------------
                ## Set statistical methods and run parameters
                ##-----------------------------------------------------------                
                max.genes=20000;max.genesets=5000
                gx.methods   = c("ttest.welch","trend.limma")
                gset.methods = c("fisher")
                extra.methods = ""
                gx.methods   = c("ttest.welch","trend.limma","edger.qlf","deseq2.wald")
                gset.methods = c("fisher","gsva","fgsea","camera","fry")
                extra.methods = c("deconv","wordcloud","connectivity")

                max.genes    = as.integer(max.genes)
                max.genesets = as.integer(max.genesets)

                ## get selected methods from input
                gx.methods    <- c(input$gene_methods,input$gene_methods2)
                gset.methods  <- c(input$gset_methods,input$gset_methods2)
                extra.methods <- c(input$extra_methods,input$extra_methods2)
                
                if(length(gx.methods)==0) {
                    shinyalert::shinyalert("ERROR","You must select at least one gene test method")
                    return(NULL)
                }
                if(length(gset.methods)==0) {
                    shinyalert::shinyalert("ERROR","You must select at least one geneset test method")
                    return(NULL)
                }

                ## at least do meta.go, infer
                extra.methods <- unique(c("meta.go","infer",extra.methods))
                
                ##----------------------------------------------------------------------
                ## Start computation
                ##----------------------------------------------------------------------
                message("[ComputePgxServer:@compute] start computations...")                
                message("[ComputePgxServer::@compute] gx.methods = ",paste(gx.methods,collapse=" "))
                message("[ComputePgxServer::@compute] gset.methods = ",paste(gset.methods,collapse=" "))
                message("[ComputePgxServer::@compute] extra.methods = ",paste(extra.methods,collapse=" "))
                
                start_time <- Sys.time()
                ## Create a Progress object
                progress <- shiny::Progress$new()
                on.exit(progress$close())    
                progress$set(message = "Processing", value = 0)
                pgx.showCartoonModal("Computation may take 5-20 minutes...")
                
                flt="";use.design=TRUE;prune.samples=FALSE
                flt <- input$filter_methods
                only.hugo <- ("only.hugo" %in% flt)
                do.protein <- ("proteingenes"   %in% flt)
                remove.unknown <- ("remove.unknown"  %in% flt)
                excl.immuno <- ("excl.immuno"  %in% flt)
                excl.xy <- ("excl.xy"  %in% flt)
                only.proteincoding <- ("only.proteincoding"  %in% flt)
                filter.genes <- ("remove.notexpressed"  %in% flt)                
               
                use.design <- !("noLM.prune" %in% input$dev_options)
                prune.samples <- ("noLM.prune" %in% input$dev_options)


                message("[ComputePgxServer:@compute] creating PGX object")                
                progress$inc(0.1, detail = "creating PGX object")            

                USE_FUTURES=1
                USE_FUTURES=0
                
                if(USE_FUTURES) {                    

                    ## !!!TRYING TO USE FUTURES. BUT SEEMS STILL TO BLOCK
                    ## OTHER SESSIONS!!!!
                    ##
                    ## IK 10.11.2021

                    message("[ComputePgxServer:@compute] using futures ")    
                    f <- future::future({
                        pgx.createPGX(
                            counts, samples, contrasts, ## genes,
                            X = NULL,   ## should we pass the pre-normalized expresson X ????
                            batch.correct = FALSE, ## done in UI                        
                            prune.samples = TRUE,  ## always prune
                            filter.genes = filter.genes,
                            ##only.chrom = FALSE,
                            ##rik.orf = !excl.rikorf,
                            only.known = !remove.unknown,
                            only.proteincoding = only.proteincoding, 
                            only.hugo = only.hugo,
                            convert.hugo = only.hugo,
                            do.cluster = TRUE,
                            cluster.contrasts = FALSE
                        )
                    })
                    ## wait until done...
                    while (!future::resolved(f)) {
                        ##cat(count, "\n")
                        message(".",appendLF = FALSE)    
                        Sys.sleep(15)  ## every 15s
                    }
                    message("done!\n")                        
                    ##value(f)
                    ngs <- future::value(f)
                    names(ngs)
                    
                } else {

                    ngs <- pgx.createPGX(
                        counts, samples, contrasts, ## genes,
                        X = NULL,   ## should we pass the pre-normalized expresson X ????
                        batch.correct = FALSE, ## done in UI                        
                        prune.samples = TRUE,  ## always prune
                        filter.genes = filter.genes,
                        ##only.chrom = FALSE,
                        ##rik.orf = !excl.rikorf,
                        only.known = !remove.unknown,
                        only.proteincoding = only.proteincoding, 
                        only.hugo = only.hugo,
                        convert.hugo = only.hugo,
                        do.cluster = TRUE,
                        cluster.contrasts = FALSE
                    )
                }
                
                names(ngs)                
                message("[ComputePgxServer:@compute] computing PGX object")
                progress$inc(0.2, detail = "computing PGX object")            

                if(USE_FUTURES) {                    
                    
                    message("[ComputePgxServer:@compute] using futures ")    
                    f <- future::future({
                        pgx.computePGX(
                            ngs,
                            max.genes = max.genes,
                            max.genesets = max.genesets, 
                            gx.methods = gx.methods,
                            gset.methods = gset.methods,
                            extra.methods = extra.methods,
                            use.design = use.design,        ## no.design+prune are combined 
                            prune.samples = prune.samples,  ##
                            do.cluster = TRUE,                
                            progress = progress,
                            lib.dir = FILES 
                        )
                    })
                    ## wait until done...
                    while (!future::resolved(f)) {
                        ##cat(count, "\n")
                        message(".",appendLF = FALSE)    
                        Sys.sleep(15)  ## every 15s
                    }
                    message("done!\n")                        
                    ##value(f)
                    ngs <- future::value(f)
                    names(ngs)

                } else {
                    
                    ngs <- pgx.computePGX(
                        ngs,
                        max.genes = max.genes,
                        max.genesets = max.genesets, 
                        gx.methods = gx.methods,
                        gset.methods = gset.methods,
                        extra.methods = extra.methods,
                        use.design = use.design,        ## no.design+prune are combined 
                        prune.samples = prune.samples,  ##
                        do.cluster = TRUE,                
                        progress = progress,
                        lib.dir = FILES 
                    )
                }
                
                end_time <- Sys.time()
                run_time  = end_time - start_time
                run_time

                message("[ComputePgxServer:@compute] total processing time of ",run_time," secs")
                
                ##----------------------------------------------------------------------
                ## annotate object
                ##----------------------------------------------------------------------
                names(ngs)
                head(ngs$samples)
                ngs$datatype = ""
                ##ngs$datatype = input$upload_datatype
                ##ngs$name = "(uploaded)"
                ngs$name = gsub("[ ]","_",input$upload_name)
                ngs$datatype    = input$upload_datatype
                ngs$description = input$upload_description
                ngs$creator <- "user"                
                
                this.date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
                ##ngs$date = date()
                ngs$date = this.date

                message("[ComputePgxServer:@compute] initialize object")
                
                ## initialize and update global PGX object
                ## ngs <- pgx.initialize(ngs)  ## here or later???
                ##uploaded$pgx <- ngs
                computedPGX(ngs)                

                ##----------------------------------------------------------------------
                ## Remove modal and show we are ready
                ##----------------------------------------------------------------------
                ##removeModal()
                if(alertready) {
                    ##beepr::beep(sample(c(3,4,5,6,8),1))  ## music!!
                    beepr::beep(2)  ## short beep
                    shinyalert::shinyalert("Ready!","We wish you lots of discoveries!")
                }

                ##for(s in names(uploaded)) uploaded[[s]] <- NULL
                ##uploaded[["pgx"]] <- NULL
                message("[ComputePgxServer:@compute] finished")

            })

            dbg("ComputePgxServer done")

            return(computedPGX)  ## pointing to reactive
        } ## end-of-server
    )
}

