##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##=====================================================================================
##============================= GADGET UI =============================================
##=====================================================================================

if(0) {

    OPG = "~/Playground/omicsplayground"
    RDIR = file.path(OPG,"R")
    FILES = file.path(OPG,"lib")
    FILESX = file.path(OPG,"libx")
    PGX.DIR = file.path(OPG,"data")    

    source("~/Playground/omicsplayground/R/pgx-include.R")
    load("~/Playground/omicsplayground/data/GSE10846-dlbcl-nc.pgx")
    source("~/Playground/omicsplayground/R/pgx-ui.R")
    source("~/Playground/omicsplayground/R/pgx-contrasts.R")    

    contr_matrix <- ngs$model.parameters$contr.matrix
    out <- gadgetize(
        makeContrastUI, makeContrastServer,
        title = "makeContrastGadget", height = 640, ## size="l", 
        pheno = ngs$samples, contr_matrix = contr_matrix )
    names(out)

    ct <- pgx.getContrasts(ngs)
    viz.Contrasts(ngs, contrasts=ct, type="volcano")

    ## Example of a shiny app
    app <- system.file("shiny-examples/bucket_list/app.R", package = "sortable")
    shiny::runApp(app)

}

MakeContrastGadget <- function(X, pheno, height=720) {
    gadgetize(MakeContrastUI, MakeContrastServer,
              title="MakeContrastGadget",
              pheno=pheno, height=height)
}

MakeContrastUI <- function(id) {
    ns <- NS(id)
    uiOutput(ns("UI"))
}

MakeContrastServerRT <- function(id, phenoRT, contrRT, countsRT=NULL, height=720)
{
    require(DT)
    moduleServer(
        id,
        function(input, output, session) {            

            dbg("*** MakeContrastServer ***\n")
            ns <- session$ns
            rv <- reactiveValues(contr=NULL)
            
            observe({
                dbg("[MakeContrastServer] obs1 : reacted\n")
                if(is.null(phenoRT()) || is.null(contrRT())) {
                    rv$contr <- NULL
                }
                req(phenoRT(), contrRT())
                rv$contr <- pgx.expMatrix(phenoRT(), contrRT())
                
                d1 <- dim(isolate(rv$contr))
                dbg("[MakeContrastServer] dim(rv)=",d1,"\n")
            })
            
            output$UI <- renderUI({
                ns <- session$ns                
                
                phenotypes <- sort(unique(c(colnames(phenoRT()),"<samples>")))
                psel <- c(grep("sample|patient|name|id|^[.]",phenotypes,value=TRUE,invert=TRUE),
                          phenotypes)[1]
                psel
                help_text = "Create comparisons by dragging conditions into the group boxes on the right."
                inline.div <- function(a) {
                    div(style="display: inline-block;vertical-align:top; width: 150px;",a)
                }
                fillCol(
                    height = 750,
                    flex = c(1,NA,NA,1),
                    fillRow(
                        flex = c(3,0.06,1.0),
                        fillCol(
                            flex = c(NA,NA,1.0),
                            h4("Create comparisons"),                        
                            ##p(help_text),
                            fillRow(
                                flex = c(1,4),
                                fillCol(
                                    flex = c(NA,NA,NA,NA,1),
                                    tipify2(
                                        selectInput(ns("param"), "Phenotype:", choices=phenotypes,
                                                    selected=psel, multiple=TRUE),
                                        "Select phenotype(s) to create conditions for making your groups."
                                    ),
                                    tipify2(
                                        textInput(ns("newname"), "Comparison name:",
                                                  placeholder="e.g. treated_vs_control"),
                                        "Give a name for your contrast. You must include _vs_ in the name to separate the names of the two groups. Try to keep the name short otherwise your plots may get ugly."),
                                    br(),
                                    tipify2(
                                        actionButton(ns("addcontrast"),"add comparison", icon=icon("plus")),
                                        "After creating the groups, press this button to add the comparison to the table."
                                    ),
                                    br()
                                ),
                                tipify(
                                    uiOutput(ns("createcomparison"),
                                             style="font-size:13px; height: 280px; overflow-y: scroll;"),
                                    "Create comparisons by dragging conditions into the main or control groups on the right. Then press add comparison to add the contrast to the table.",
                                    placement="top", options = list(container = "body"))
                            )
                        ),
                        br(),
                        ##plotOutput(ns("pcaplot"), height="330px")
                        plotWidget(ns("pcaplot"))
                    ),
                    h4("Contrast table"),
                    fillRow(
                        height = 24,
                        flex = c(NA,1),
                        actionButton(ns("autocontrast"),"add auto-contrasts", icon=icon("plus"),
                                     class="small-button"),
                        br()
                    ),
                    ##tags$head(tags$style("table.dataTable.compact tbody th, table.dataTable.compact tbody td {padding: 0px 10px;}")),
                    ## this.style(ns("contrastTable"), "table.dataTable.compact tbody th, table.dataTable.compact tbody td {padding: 0px 10px;}"),                    
                    div(DT::dataTableOutput(ns("contrastTable")),
                        style="font-size:13px; height: 300px; margin-top: 10px;overflow-y: scroll;")
                )
                
            })
            outputOptions(output, "UI", suspendWhenHidden=FALSE) ## important!!!
            
            sel.conditions <- reactive({
                req(phenoRT())
                df <- phenoRT()
                df$"<samples>" <- rownames(df)
                pp <- intersect(input$param, colnames(df))
                cond <- apply(df[,pp,drop=FALSE],1,paste,collapse="_")
                cond
            })
            
            output$createcomparison <- renderUI({

                library(sortable)
                req(input$param)
                cond <- sel.conditions()
                if(length(cond)==0) return(NULL)
                items <- c("<others>",sort(unique(cond)))
                dbg("[MakeContrastServer:createcomparison] items=",items,"\n")
                
                tagList(
                    tags$head(tags$style(".default-sortable .rank-list-item {padding: 2px 15px;}")),
                    bucket_list(
                        ##header = h4("Create comparison:"),
                        header = NULL,
                        add_rank_list(
                            text = "Conditions:",
                            labels = items
                        ),
                        add_rank_list(
                            input_id = ns("group1"),
                            text = "Main group:"
                        ),
                        add_rank_list(
                            input_id = ns("group2"),
                            text = "Control group:"
                        ),
                        group_name = "cmpbucket"
                    )
                )
            })
            
            buttonInput <- function(FUN, len, id, ...) {
                inputs <- character(len)
                for (i in seq_len(len)) {
                    inputs[i] <- as.character(FUN(paste0(id, i), ...))
                }
                inputs
            }    

            observeEvent( c(input$group1, input$group2), {
                g1 <- gsub("[-_.,<> ]",".",input$group1)
                g2 <- gsub("[-_.,<> ]",".",input$group2)
                g1 <- gsub("[.]+",".",g1)
                g2 <- gsub("[.]+",".",g2)                
                g1 <- paste(g1,collapse="")
                g2 <- paste(g2,collapse="")
                if(is.null(g1) || length(g1)==0) g1 <- ""
                if(is.null(g2) || length(g2)==0) g2 <- ""
                if(is.na(g1)) g1 <- ""
                if(is.na(g2)) g2 <- ""
                g1 <- substring(g1,1,20)
                g2 <- substring(g2,1,20)
                pp <- paste(input$param,collapse=".")
                pp <- gsub("[-_.,<> ]","",pp)
                tt <- paste0(pp,":",g1,"_vs_",g2)
                if(g1=="" && g2=="") tt <- ""
                updateTextInput(session, "newname", value=tt)
            })
            
            observeEvent( input$contrast_delete, {
                ## Observe if a contrast is to be deleted
                ##
                id <- as.numeric(gsub(".*_","",input$contrast_delete))
                message('[contrast_delete] clicked on delete contrast',id)
                if(length(id)==0) return(NULL)
                ##updateActionButton(session, paste0("contrast_delete_",id),label="XXX")
                message("[contrast_delete] 1: dim(ct) = ",dim(rv$contr))
                if(!is.null(rv$contr) && NCOL(rv$contr) <= 1) {
                    rv$contr <- rv$contr[,0,drop=FALSE]
                } else {
                    rv$contr <- rv$contr[,-id,drop=FALSE] 
                }
                message("[contrast_delete] 2: dim(ct) = ",dim(rv$contr))
                ## invalidateLater(3000, session)
            })
            
            observeEvent( input$addcontrast, {
                
                cond <- sel.conditions()
                in.main <- cond %in% input$group1
                in.ref1 <- cond %in% input$group2 
                in.ref2 <- ("<others>" %in% input$group2) &
                    (!cond %in% input$group1)
                in.ref <- in.ref1 | in.ref2                
                ctx <- 1*(in.main) - 1*(in.ref)

                ##ct.name <- paste0(input$group1name,"_vs_",input$group2name)
                ct.name <- input$newname
                samples <- rownames(phenoRT())

                if( sum(ctx<0)==0 || sum(ctx>0)==0 ) {
                    shinyalert("ERROR","Both groups must have samples")
                    return(NULL)
                }
                if(ct.name %in% c(NA,""," ")) {
                    shinyalert("ERROR","You must give a contrast name")
                    return(NULL)
                }
                if(!is.null(rv$contr) && ct.name %in% colnames(rv$contr)) {
                    shinyalert("ERROR","Contrast name already exists. Please give a different name.")
                    return(NULL)
                }
                
                ## update reactive value
                if(!is.null(rv$contr)) {
                    rv$contr <- cbind(rv$contr, ctx)
                    colnames(rv$contr)[ncol(rv$contr)] <- ct.name
                } else {
                    rv$contr <- matrix(ctx, ncol=1, dimnames=list(samples,ct.name))
                }
                
            })

            observeEvent( input$autocontrast, {

                req(phenoRT())
                df <- phenoRT()
                dbg("[MakeContrastServerRT:autocontrast] dim(df)=",dim(df),"\n")                
                ct <- pgx.makeAutoContrast(
                    df, mingrp=3, slen=20, ref=NULL, fix.degenerate=FALSE)
                rownames(ct$exp.matrix) <- rownames(df)
                ctx <- ct$exp.matrix
                
                dbg("[MakeContrastServerRT:autocontrast] updating contrasts...\n")
                dbg("[MakeContrastServerRT:autocontrast] dim(rv$contr)=",dim(rv$contr),"\n")
                dbg("[MakeContrastServerRT:autocontrast] dim(ctx)=",dim(ctx),"\n")

                ## update reactive value
                if(!is.null(rv$contr)) {
                    rv$contr <- cbind(rv$contr, ctx)
                } else {
                    rv$contr <- ctx
                }

            })

            output$contrastTable <- DT::renderDataTable({

                message('[contrastTable] >>>> reacted')
                message('[contrastTable] dim(contr) = ',dim(rv$contr))
                
                ct <- rv$contr
                
                if(is.null(ct) || NCOL(ct)==0) {
                    df <- data.frame(
                        delete = 0,
                        comparison = "",
                        n1 = 0,
                        n2 = 0,
                        samples1 = "",
                        samples2 = ""
                    )[0,]
                } else {
                    paste.max <- function(x,n=6) {
                        ##x <- unlist(x)
                        if(length(x)>n) {
                            x <- c(x[1:n], paste("+",length(x)-n,"others"))
                        }
                        paste(x,collapse=" ")
                    }                    
                    ct1  <- sign(ct)
                    if(NCOL(ct)==1) {
                        ss1 <- names(which(ct1[,1]==1))
                        ss2 <- names(which(ct1[,1]==-1))
                        ss1 <- paste.max(ss1,6)
                        ss2 <- paste.max(ss2,6)
                    } else {
                        ss0 <- rownames(ct)
                        ss1 <- apply(ct1,2,function(x) paste.max(ss0[which(x==1)])) 
                        ss2 <- apply(ct1,2,function(x) paste.max(ss0[which(x==-1)]))
                    }

                    deleteButtons <- buttonInput(
                        FUN = actionButton,
                        len = ncol(ct),
                        ##id = 'contrast_delete_',
                        id = paste0('contrast_delete_',sample(99999,1),"_"),  ## hack to allow double click
                        label = "delete",
                        ##size = "mini",
                        width = "50px",
                        inline = TRUE,
                        icon = icon("trash-alt"),
                        class = "btn-inline",
                        style='padding:2px; margin:2px; font-size:95%; color: #B22222;',
                        ##onclick = 'Shiny.onInputChange(\"contrast_delete\",this.id)'
                        onclick = paste0('Shiny.onInputChange(\"',ns("contrast_delete"),'\",this.id)')
                    )

                    df <- data.frame(
                        delete = deleteButtons,
                        comparison = colnames(ct),
                        n1 = colSums(ct > 0),
                        n2 = colSums(ct < 0 ),
                        samples1 = ss1,
                        samples2 = ss2
                    )
                }
                rownames(df) <- NULL
                dbg("[contrastTable] render datatable dim(df) = ",dim(df))
                
                datatable(
                    df, rownames=FALSE,
                    escape = c(-1),
                    selection = 'none',
                    class="compact cell-border",
                    options = list(
                        dom = 't',
                        pageLength = 999,
                        ## autoWidth = TRUE, ## scrollX=TRUE,
                        columnDefs = list(
                            list(width='20px', targets=c(0,2,3)),
                            list(width='150px', targets=c(1)),
                            list(width='400px', targets=c(4,5))
                        )                        
                    )
                ) %>%
                    DT::formatStyle(0, target='row', fontSize='12px', lineHeight='99%')                
            }, server=FALSE)
            

            pcaplot.RENDER <- reactive({

                ##ngs <- inputData()
                ##X <- ngs$X
                pheno <- phenoRT()
                counts <- countsRT()
                if(is.null(pheno) || is.null(counts)) return(NULL)
                if(NCOL(pheno)==0 || NCOL(counts)==0) return(NULL)
                req(pheno)
                req(counts)
                
                method <- input$pcaplot.method
                X <- log2(1 + counts)
                clust <- pgx.clusterMatrix(X, dims=2, method=method)
                names(clust)

                y <- sel.conditions()
                
                par(mar=c(4,1,1,1))
                pgx.scatterPlotXY(
                    clust$pos2d, var=y, plotlib="plotly",
                    legend = FALSE ##, labels=TRUE
                )
                
            })
            
            pcaplot.opts = tagList(
                tipify( selectInput( ns("pcaplot.method"), "Method:",
                                    choices = c("pca","tsne","umap"),
                                    width = '100%'),"Choose clustering method.",
                       placement="right", options = list(container = "body"))
            )

            observe({
                req(phenoRT())
                px <- colnames(phenoRT())
                message("[MakeContrastModule] px = ", paste(px,collapse=" "))
                updateSelectInput(session, "pcaplot.colvar", choices=px)
            })
            
            callModule(
                plotModule, 
                id = "pcaplot",
                func = pcaplot.RENDER, ## ns=ns,
                plotlib = "plotly", 
                options = pcaplot.opts,
                height = c(320,700), width=c("auto",800),
                pdf.width=8, pdf.height=8,
                title="PCA/tSNE plot",
                ##info.text = hm_PCAplot_text
                ##caption = pca_caption_static
            )

            
            ##ct <- rv$contr
            return(reactive({rv$contr}))  ## pointing to reactive
        } ## end-of-server
    )
    
}

