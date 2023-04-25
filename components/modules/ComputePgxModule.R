##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

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

ComputePgxServer <- function(
    id,
    countsRT,
    samplesRT,
    contrastsRT,
    batchRT,
    metaRT,
    lib.dir,
    pgx.dirRT,
    enable_button = TRUE,
    alertready = TRUE,
    max.genes = 20000,
    max.genesets = 10000,
    max.datasets = 100,
    height = 720,
    r_global)
{
    shiny::moduleServer(
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
            EXTRA.METHODS = c("deconv", "drugs", "wordcloud","connectivity", "wgcna")
            EXTRA.NAMES = c("celltype deconvolution", "drugs connectivity",
                            "wordcloud","experiment similarity", "WGCNA")
            EXTRA.SELECTED = c("deconv","drugs","wordcloud","connectivity","wgcna")

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
                                                    class="btn-outline-primary"),
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
                                style = "width: 95%;",
                                shiny::checkboxGroupInput(
                                    ns('filter_methods'),
                                    shiny::HTML('<h4>Feature filtering:</h4><br/>'),
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
                                style = "width: 95%;",
                                shiny::checkboxGroupInput(
                                    ns('gene_methods'),
                                    shiny::HTML('<h4>Gene tests:</h4><br/>'),
                                    GENETEST.METHODS,
                                    selected = GENETEST.SELECTED
                                )
                            ),
                            shiny::wellPanel(
                                style = "width: 95%;",
                                shiny::checkboxGroupInput(
                                    ns('gset_methods'),
                                    shiny::HTML('<h4>Enrichment methods:</h4><br/>'),
                                    GENESET.METHODS,
                                    selected = GENESET.SELECTED
                                ),
                            ),
                            shiny::wellPanel(
                                style = "width: 95%;",
                                shiny::checkboxGroupInput(
                                    ns('extra_methods'),
                                    shiny::HTML('<h4>Extra analysis:</h4><br/>'),
                                    choiceValues = EXTRA.METHODS,
                                    choiceNames = EXTRA.NAMES,
                                    selected = EXTRA.SELECTED
                                )
                            ),
                            shiny::wellPanel(
                                shiny::checkboxGroupInput(
                                    ns('dev_options'),
                                    shiny::HTML('<h4>Developer options:</h4><br/>'),
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

                meta <- metaRT()

                if(!is.null(meta[['name']])) {
                    shiny::updateTextInput(session, "upload_name", value=meta[['name']])
                    ## shiny::updateTextInput(session, ns("upload_name"), value=meta[['name']])
                }
                if(!is.null(meta[['description']])) {
                    shiny::updateTextAreaInput(session, "upload_description", value=meta[['description']])
                    ##shiny::updateTextAreaInput(session, ns("upload_description"), value=meta[['description']])
                }

            })

            ##------------------------------------------------------------------
            ## After confirmation is received, start computing the PGX
            ## object from the uploaded files
            ## ------------------------------------------------------------------

            # Define a reactive value to store the process object
            process_obj  <- reactiveVal(NULL)
            computedPGX  <- shiny::reactiveVal(NULL)
            temp_dir     <- reactiveVal(NULL)
            process_counter <- reactiveVal(0)
            reactive_timer  <- reactiveTimer(20000)  # Triggers every 10000 milliseconds (20 second)

            shiny::observeEvent( input$compute, {
                ## shiny::req(input$upload_hugo,input$upload_filtergenes)

                ##-----------------------------------------------------------
                ## Check validity
                ##-----------------------------------------------------------
                if(!enable_button()) {
                    message("[ComputePgxServer:@compute] WARNING:: *** NOT ENABLED ***")
                    return(NULL)
                }

                pgxdir <- pgx.dirRT()

                numpgx <- length(dir(pgxdir, pattern="*.pgx$"))

                if(numpgx >= max.datasets) {
                    ### should use sprintf here...
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
                extra.methods = c("deconv","wordcloud","connectivity", "wgcna")

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

                # start_time <- Sys.time()
                ## Create a Progress object
                # progress <- shiny::Progress$new()
                # progress$set(message = "Processing", value = 0)
                # pgx.showCartoonModal("Computation may take 5-20 minutes...")

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

                # create folder with random name to store the csv files

                # Generate random name for temporary folder

                # Create temporary folder
                temp_dir(tempfile(pattern = "pgx_"))
                dir.create(temp_dir())
                dbg("[compute PGX process] : tempFile", temp_dir())

                this.date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
                path_to_params <- file.path(temp_dir(), "params.RData")
                dataset_name <- gsub("[ ]","_",input$upload_name)

                # Define create_pgx function arguments
                params <- list(
                    samples = samples,
                    counts = counts,
                    contrasts = contrasts,
                    batch.correct = FALSE,
                    prune.samples = TRUE,
                    filter.genes = filter.genes,
                    only.known = !remove.unknown,
                    only.proteincoding = only.proteincoding,
                    only.hugo = only.hugo,
                    convert.hugo = only.hugo,
                    do.cluster = TRUE,
                    cluster.contrasts = FALSE,
                    max.genes = max.genes,
                    max.genesets = max.genesets,
                    gx.methods = gx.methods,
                    gset.methods = gset.methods,
                    extra.methods = extra.methods,
                    use.design = use.design,        ## no.design+prune are combined
                    prune.samples = prune.samples,  ##
                    do.cluster = TRUE,
                    libx.dir = lib.dir, # needs to be replaced with libx.dir
                    name = dataset_name,
                    datatype = input$upload_datatype,
                    description = input$upload_description,
                    creator = "user",
                    this.date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    date = this.date
                )

                ## override username from session object
                params$creator <- session$user


                saveRDS(params, file=path_to_params)

                # Normalize paths
                script_path <- normalizePath(file.path(get_opg_root(), "bin", "pgxcreate_op.R"))
                cmd <- normalizePath(temp_dir())

                # Start the process and store it in the reactive value
                shinyalert::shinyalert(
                    title = "Crunching your data!",
                    text = paste0("Your dataset will be computed in the background. You can continue to play with a different dataset in the meantime. When it is ready, it will appear in your dataset library.")
                    ## timer = 8000
                )
                bigdash.selectTab(
                    session,
                    selected = 'load-tab'
                )

                dbg("[compute PGX process] : starting process")
                process_counter(process_counter() + 1)

                if(is.null(process_obj())) {
                    process_obj(
                        list(
                            list(
                                process = processx::process$new(
                                     "Rscript",
                                     args = c(script_path, cmd),
                                     supervise = TRUE,
                                     stderr = '|',
                                     stdout = '|'
                                ),
                                dataset_name = gsub("[ ]","_",input$upload_name),
                                temp_dir = temp_dir())))
                } else {
                    process_obj(
                        append(
                            process_obj(),
                            list(
                                list(
                                    process = processx::process$new(
                                        "Rscript",
                                        args = c(script_path, cmd),
                                        supervise = TRUE,
                                        stderr = '|',
                                        stdout = '|'
                                    ),
                                    dataset_name = gsub("[ ]","_",input$upload_name),
                                    temp_dir = temp_dir())
                                ))
                            )
                }
            })

            check_process_status <- reactive({
                if (process_counter() == 0) {
                    return(NULL)
                }

                reactive_timer()
                active_processes <- process_obj()
                completed_indices <- c()

                for (i in seq_along(active_processes)) {
                    #i=1
                    current_process <- active_processes[[i]]$process
                    temp_dir <- active_processes[[i]]$temp_dir

                    process_status <- current_process$get_exit_status()
                    process_alive <- current_process$is_alive()

                    if (!is.null(process_status) && process_status == 0) {
                        # Process completed successfully
                        dbg("[compute PGX process] : process completed")
                        on_process_completed(temp_dir = temp_dir)
                        completed_indices <- c(completed_indices, i)
                    } else if (!is.null(process_status) && process_status != 0) {
                        on_process_error()
                        completed_indices <- c(completed_indices, i)
                    } else {
                        # Process is still running, do nothing
                        dbg("[compute PGX process] : process still running")

                        ## write error to console and temp file
                        stderr_output <- current_process$read_error_lines()
                        ##logfile <- file.path(temp_dir,"process.log")
                        logfile <- normalizePath(file.path(OPG,"processx.log"))
                        dbg("[compute PGX process] : writing stderr to ", logfile)
                        writeLines(stderr_output, logfile)
                    }

                }

                # Remove completed processes from the list
                if (length(completed_indices) > 0) {
                    active_processes <- active_processes[-completed_indices]
                    process_obj(active_processes)
                }

                return(NULL)
            })

            # Function to execute when the process is completed successfully
            on_process_completed <- function(temp_dir) {
                dbg("[compute PGX process] on_process_completed() called!")
                process_counter(process_counter()-1) # stop the timer
                result_pgx <- file.path(temp_dir, "my.pgx")
                if (file.exists(result_pgx)) {
                    load(result_pgx)  ## always pgx
                    computedPGX(pgx)
                } else {
                    message("[compute PGX process] : Error: Result file not found")
                }
                unlink(temp_dir, recursive = TRUE)
            }

            on_process_error <- function() {
                dbg("[compute PGX process] on_process_error() called!")
                process_counter(process_counter()-1) # stop the timer
                message("Error: Process completed with an error")

                proc <-process_obj()[[1]]$process
                stderr_output <- proc$read_error_lines()

                if (length(stderr_output) > 0) {
                    message("Standard error output from the process:")
                    ##for (line in stderr_output) { message(line) }
                    writeLines(stderr_output, con=stderr())
                } else {
                    message("No standard error output available from the process")
                }
            }

            ## what does this do???
            observe(check_process_status())

            observe({
                if (process_counter() > 0){
                    shiny::insertUI(
                        selector = ".current-dataset",
                        where = "beforeEnd",
                        ui = loading_spinner("Computation in progress...")
                        )
                } else if (process_counter() == 0) {
                    shiny::removeUI(selector = ".current-dataset > #spinner-container")
                }

                if (process_counter() < opt$MAX_DS_PROCESS) {
                    shinyjs::enable("compute")
                } else {
                    shinyjs::disable("compute")
                }
                })

            return(computedPGX)
        } ## end-of-server
    )
}

