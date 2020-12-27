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
                             alertready=TRUE, max.features = 10000)
{
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns

            ## statistical method for GENE level testing
            GENETEST.METHODS = c("ttest","ttest.welch","voom.limma","trend.limma","notrend.limma",
                                 "deseq2.wald","deseq2.lrt","edger.qlf","edger.lrt")
            GENETEST.METHODS1 = c("ttest","ttest.welch","voom.limma","trend.limma","notrend.limma")
            GENETEST.METHODS1 = c("ttest","ttest.welch","voom.limma","trend.limma","notrend.limma")
            GENETEST.METHODS2 = c("deseq2.wald","deseq2.lrt","edger.qlf","edger.lrt")

            ## statistical method for GENESET level testing
            GENESET.METHODS = c("fisher","ssgsea","gsva", "spearman", "camera", "fry",
                                ##"plage","enricher","gsea.permPH","gsea.permGS","gseaPR",
                                "fgsea")
            GENESET.METHODS1 = c("fisher","spearman","fgsea") 
            GENESET.METHODS2 = c("gsva","ssgsea", "camera", "fry")

            ## batch correction and extrs methods
            BC_METHODS = c("limma","SVA","ComBat")
            EXTRA_METHODS = c("meta.go","infer","deconv","drugs-combo",
                              "wordcloud","connectivity")
            EXTRA_METHODS1 = c("meta.go","infer","deconv","drugs","wordcloud")
            EXTRA_METHODS2 = c("drugs-combo","connectivity")
            

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
                                          "excl.rikorf"
                                          ## "excl.immuno"
                                          ## "excl.xy"
                                          ),
                                    selected = c("only.hugo","only.proteincoding",
                                                 "excl.rikorf"),
                                    choiceNames =
                                        c("Convert to HUGO",
                                          "Protein-coding only",
                                          "Exclude Rik/ORF"
                                          ##"Exclude immunogenes",
                                          ##"Exclude X/Y genes"
                                          )
                                )
                            ),
                            wellPanel(
                                checkboxGroupInput(
                                    ns('gene_methods'),
                                    'Gene tests:',
                                    GENETEST.METHODS1,
                                    selected = c("ttest","ttest.welch","trend.limma")
                                ),
                                premium.feature(
                                    checkboxGroupInput(
                                        ns('gene_methods2'),
                                        NULL,
                                        GENETEST.METHODS2
                                    )
                                )
                            ),
                            wellPanel(
                                checkboxGroupInput(
                                    ns('gset_methods'),
                                    'Enrichment methods:',
                                    GENESET.METHODS1,
                                    selected = c("fisher","spearman","fgsea")
                                ),
                                premium.feature(
                                    checkboxGroupInput(
                                        ns('gset_methods2'),
                                        NULL,
                                        GENESET.METHODS2
                                    )
                                )
                            ),
                            wellPanel(
                                checkboxGroupInput(
                                    ns('extra_methods'),
                                    'Extra analysis:',
                                    EXTRA_METHODS1,
                                    selected = c("infer","wordcloud")
                                ),
                                premium.feature(
                                    checkboxGroupInput(
                                        ns('extra_methods2'),
                                        NULL,
                                        EXTRA_METHODS2
                                    )
                                )
                                
                            ),
                            div(width=150)                            
                        ) ## end of fillRow
                    ) ## end of conditional panel
                ) ## end of fill Col
            })


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
                
                cat("[LoadingModule] dim(counts)=",dim(counts),"\n")
                cat("[LoadingModule] dim(samples)=",dim(samples),"\n")
                cat("[LoadingModule] dim(contrasts)=",dim(contrasts),"\n")
                
                ##-----------------------------------------------------------
                ## Set statistical methods and run parameters
                ##-----------------------------------------------------------                
                max.genes    = as.integer(max.features)
                max.genesets = as.integer(max.features)
                
                gx.methods   = c("ttest.welch","trend.limma","edger.qlf","deseq2.wald")
                gset.methods = c("fisher","gsva","fgsea","camera","fry")
                extra.methods = c("meta.go","infer","deconv","drugs-combo",
                                  "wordcloud","connectivity")

                ## get selected methods from input
                gx.methods   <- c(input$gene_methods,input$gene_methods2)
                gset.methods <- c(input$gset_methods,input$gset_methods2)
                extra.methods <- c(input$extra_methods,input$extra_methods2)
                
                ##----------------------------------------------------------------------
                ## Start computation
                ##----------------------------------------------------------------------
                start_time <- Sys.time()
                ## Create a Progress object
                progress <- shiny::Progress$new()
                on.exit(progress$close())    
                progress$set(message = "Processing", value = 0)

                flt <- input$filter_methods
                do.filter <- (length(flt)>0 && flt[1]!="")
                only.hugo <- ("only.hugo" %in% flt)
                do.protein <- ("proteingenes"   %in% flt)
                excl.rikorf <- ("excl.rikorf"  %in% flt)
                excl.immuno <- ("excl.immuno"  %in% flt)
                excl.xy <- ("excl.xy"  %in% flt)
                only.chrom <- ("only.chrom"  %in% flt)
                only.proteincoding <- ("only.proteincoding"  %in% flt)                
                
                progress$inc(0.01, detail = "parsing data")            
                ngs <- pgx.createPGX(
                    counts, samples, contrasts, ## genes,
                    X = NULL,
                    batch.correct = FALSE,      ## done in UI                        
                    filter.genes = do.filter,
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

