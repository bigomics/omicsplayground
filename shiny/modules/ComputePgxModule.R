##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

if(0) {
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
        title="ComputePGX",
        countsRT = reactive(counts),
        samplesRT = reactive(samples),
        contrastsRT = reactive(contrasts)
    )
}

ComputePgxUI <- function(id) {
    ns <- NS(id)
    uiOutput(ns("UI"))
}

ComputePgxServer <- function(id, countsRT, samplesRT, contrastsRT,
                             FILES, enable=TRUE, height=720,
                             alertready=TRUE,
                             max.genes = 20000, max.genesets = 10000)
{
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns

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
            EXTRA.NAMES = c("celltype deconvolution", "drugs enrichment",
                            "wordcloud","experiment connectivity")
            EXTRA.SELECTED = c("drugs","wordcloud")
            
            output$UI <- renderUI({
                fillCol(
                    height = height,
                    flex = c(0.2,NA,0.05,1.3),
                    br(),
                    ## conditionalPanel(
                    ##     "input.options%2 == 0", ns=ns,
                    ##     div(style="height: 250px;")
                    ## ),
                    fluidRow(
                        column(
                            4, align="center", offset=4,
                            textInput(
                                ns("upload_name"),NULL, ##"Dataset:",
                                ##width="100%",
                                ## width=420,
                                placeholder="Name of your dataset"),
                            br(), 
                            div(textAreaInput(
                                ns("upload_description"), NULL, ## "Description:",
                                placeholder="Give a short description of your dataset",
                                ##width="100%",
                                ## width=400
                                height=100, resize='none'),
                                style="margin-left: 0px;"),
                            br(),
                            ##textInput(ns("upload_datatype"),"Datatype:"), br(),
                            actionButton(ns("compute"),"Compute!",icon=icon("running"),
                                         class="run-button"),
                            br(),br(),
                            actionLink(ns("options"), "Advanced", icon=icon("cog", lib="glyphicon")),
                            style = ""
                        )
                    ),
                    br(),
                    conditionalPanel(
                        "input.options%2 == 1", ns=ns,
                        fillRow(
                            div(width=150),
                            wellPanel(
                                checkboxGroupInput(
                                    ns('filter_methods'),
                                    'Feature filtering:',
                                    choiceValues =
                                        c("only.hugo",
                                          "only.proteincoding",
                                          "excl.rikorf",
                                          "remove.notexpressed"
                                          ## "excl.immuno"
                                          ## "excl.xy"
                                          ),
                                    choiceNames =
                                        c("convert to HUGO",
                                          "protein-coding only",
                                          "exclude Rik/ORF",
                                          "remove not-expressed"
                                          ##"Exclude immunogenes",
                                          ##"Exclude X/Y genes"
                                          ),
                                    selected = c(
                                        "only.hugo",
                                        "only.proteincoding",
                                        "excl.rikorf",
                                        "remove.notexpressed"
                                    )
                                )
                            ),
                            wellPanel(
                                checkboxGroupInput(
                                    ns('gene_methods'),
                                    'Gene tests:',
                                    GENETEST.METHODS,
                                    selected = GENETEST.SELECTED
                                )
                            ),
                            wellPanel(
                                checkboxGroupInput(
                                    ns('gset_methods'),
                                    'Enrichment methods:',
                                    GENESET.METHODS,
                                    selected = GENESET.SELECTED
                                ),
                            ),
                            wellPanel(
                                checkboxGroupInput(
                                    ns('extra_methods'),
                                    'Extra analysis:',
                                    choiceValues = EXTRA.METHODS,
                                    choiceNames = EXTRA.NAMES,
                                    selected = EXTRA.SELECTED
                                )
                            ),
                            div(width=150)                            
                        ) ## end of fillRow
                    ) ## end of conditional panel
                ) ## end of fill Col
            })

            if(FALSE) {
                observeEvent( input$gene_methods, {
                    if(length(input$gene_methods) > 3){
                        updateCheckboxGroupInput(session, "gene_methods",
                                                 selected= tail(input$gene_methods,3))
                    }
                    if(length(input$gene_methods) < 1){
                        updateCheckboxGroupInput(session, "gene_methods", selected= "ttest")
                    }
                })
            }

            
            observeEvent( input$options, {
                ## shinyjs::disable(ns("gene_methods2"))
            })
            
            observeEvent(enable(), {
                ## NEED CHECK. not working... 
                ##
                cat("[ComputePgxServer:@enable] reacted\n")
                if(!enable()){
                    cat("[ComputePgxServer:@enable] disabling compute button\n")
                    shinyjs::disable(ns("compute"))
                } else {
                    cat("[ComputePgxServer:@enable] enabling compute button\n")
                    shinyjs::enable(ns("compute"))
                }
            })

            output$correction_summary <- renderText({                
            })
                        
            ##------------------------------------------------------------------
            ## After confirmation is received, start computing the PGX
            ## object from the uploaded files
            ## ------------------------------------------------------------------
            computedPGX  <- reactiveVal(NULL)
            
            observeEvent( input$compute, {
                
                cat("[ComputePgxServer::@compute] reacted!\n")
                ## req(input$upload_hugo,input$upload_filtergenes)
                
                if(!enable()) {
                    cat("[ComputePgxServer:@compute] *** NOT ENABLED ***\n")
                    return(NULL)
                }

                has.name <- input$upload_name != ""
                has.description <- input$upload_description != ""
                if(!has.name || !has.description) {
                    shinyalert("ERROR","You must give a dataset name and description")
                    return(NULL)
                }
                
                ##-----------------------------------------------------------
                ## Retrieve the most recent matrices from reactive values
                ##-----------------------------------------------------------                
                counts    <- countsRT()
                samples   <- samplesRT()
                samples   <- data.frame(samples, stringsAsFactors=FALSE,check.names=FALSE)
                contrasts <- as.matrix(contrastsRT())
                contrasts[is.na(contrasts)] <- 0
                
                dbg("[LoadingModule] dim(counts)=",dim(counts))
                dbg("[LoadingModule] dim(samples)=",dim(samples))
                dbg("[LoadingModule] dim(contrasts)=",dim(contrasts))
                
                ##-----------------------------------------------------------
                ## Set statistical methods and run parameters
                ##-----------------------------------------------------------                
                max.genes    = as.integer(max.genes)
                max.genesets = as.integer(max.genesets)
                
                gx.methods   = c("ttest.welch","trend.limma","edger.qlf","deseq2.wald")
                gset.methods = c("fisher","gsva","fgsea","camera","fry")
                extra.methods = c("deconv","wordcloud","connectivity")

                ## get selected methods from input
                gx.methods   <- c(input$gene_methods,input$gene_methods2)
                gset.methods <- c(input$gset_methods,input$gset_methods2)
                extra.methods <- c(input$extra_methods,input$extra_methods2)

                if(length(gx.methods)==0) {
                    shinyalert("ERROR","You must select at least one gene test method")
                    return(NULL)
                }
                if(length(gset.methods)==0) {
                    shinyalert("ERROR","You must select at least one geneset test method")
                    return(NULL)
                }

                ## at least do meta.go, infer
                extra.methods <- unique(c("meta.go","infer",extra.methods))
                
                ##----------------------------------------------------------------------
                ## Start computation
                ##----------------------------------------------------------------------
                start_time <- Sys.time()
                ## Create a Progress object
                progress <- shiny::Progress$new()
                on.exit(progress$close())    
                progress$set(message = "Processing", value = 0)

                pgx.showCartoonModal("Computation may take 5-20 minutes...")
                
                flt <- input$filter_methods
                only.hugo <- ("only.hugo" %in% flt)
                do.protein <- ("proteingenes"   %in% flt)
                excl.rikorf <- ("excl.rikorf"  %in% flt)
                excl.immuno <- ("excl.immuno"  %in% flt)
                excl.xy <- ("excl.xy"  %in% flt)
                only.chrom <- ("only.chrom"  %in% flt)
                only.proteincoding <- ("only.proteincoding"  %in% flt)
                filter.genes <- ("remove.notexpressed"  %in% flt)                
                
                progress$inc(0.01, detail = "parsing data")            
                ngs <- pgx.createPGX(
                    counts, samples, contrasts, ## genes,
                    X = NULL,
                    batch.correct = FALSE,      ## done in UI                        
                    filter.genes = filter.genes,
                    only.chrom = only.chrom,
                    rik.orf = !excl.rikorf,
                    only.proteincoding = only.proteincoding, 
                    only.hugo = only.hugo,
                    convert.hugo = only.hugo
                )
                names(ngs)

                if(excl.xy) {
                    ## fill me
                }

                ngs <- pgx.computePGX(
                    ngs,
                    max.genes = max.genes,
                    max.genesets = max.genesets, 
                    gx.methods = gx.methods,
                    gset.methods = gset.methods,
                    extra.methods = extra.methods,
                    lib.dir = FILES, do.cluster=TRUE,                
                    progress = progress)
                
                end_time <- Sys.time()
                run_time  = end_time - start_time
                run_time

                message("upload_compute :: total processing time of ",run_time," secs")
                
                ##----------------------------------------------------------------------
                ## annotate object
                ##----------------------------------------------------------------------
                names(ngs)
                head(ngs$samples)
                ngs$datatype = ""
                ##ngs$datatype = input$upload_datatype
                ##ngs$name = "(uploaded)"
                ngs$name = gsub("[ ]","_",input$upload_name)
                ngs$description = input$upload_description

                this.date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
                ##ngs$date = date()
                ngs$date = this.date
                
                ## initialize and update global PGX object
                ngs <- pgx.initialize(ngs)
                ##uploaded$pgx <- ngs
                computedPGX(ngs)                

                ##----------------------------------------------------------------------
                ## Remove modal and show we are ready
                ##----------------------------------------------------------------------
                ##removeModal()
                if(alertready) {
                    shinyalert("Ready!","We wish you lots of discoveries!")
                }

                ##for(s in names(uploaded)) uploaded[[s]] <- NULL
                ##uploaded[["pgx"]] <- NULL

            })

            return(computedPGX)  ## pointing to reactive
        } ## end-of-server
    )
}

