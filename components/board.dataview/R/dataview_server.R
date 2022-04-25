##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' DataView module server function
#'
#' @description A shiny Module (server code).
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#' @param pgx Reactive expression that provides the input pgx data object 
#'
#' @export 
DataViewBoard <- function(id, pgx)
{
  moduleServer(id, function(input, output, session)
  {
    ns <- session$ns ## NAMESPACE
    rowH = 355  ## row height of panels
    imgH = 315  ## height of images
    fullH = 750 ## full height of panel
    tabH = 600  ## height of tables

    ##----------------------------------------------------------------------
    ## More Info (pop up window)
    ##----------------------------------------------------------------------
    dropdown_search_gene='<code>Search gene</code>'
    menu_grouped='<code>Group by</code>'

    data_infotext =paste0(
        'The <strong>DataView module</strong> provides information and visualisations of the dataset to quickly lookup a gene,
        check the counts, or view the data tables.<br><br>
        The <strong>Plots</strong> panel displays figures related to the expression level of the selected gene,
        correlation, and average expression ranking within the dataset.
        More information about the gene and hyperlinks to external databases are provided. Furthermore,
        it displays the correlation and tissue expression for a selected gene in external reference datasets.
        In the <strong>Counts</strong> panel, the total number of counts (abundance) per sample and their distribution among the samples are displayed.
        This is most useful to check the technical quality of the dataset, such as total read counts or abundance of ribosomal genes.
        In <strong>Gene Table</strong> panel, the exact expression values across the samples can be looked up,
        where genes are ordered by the correlation with respect to the first gene. Gene-wise average expression of a phenotype sample grouping
        is also presented in this table. In the <strong>Samples</strong> panel, more complete information about samples can be found.
        Finally, the <strong>Contrasts</strong> panel, shows information about the phenotype comparisons.
        <br><br><br>
        <center><iframe width="560" height="315" src="https://www.youtube.com/embed/S32SPINqO8E"
        title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture"
        allowfullscreen></iframe></center>
    ')

   
    ## ------- observe functions -----------
    shiny::observeEvent( input$data_info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>Data View Board</strong>"),
            shiny::HTML(data_infotext),
            easyClose = TRUE, size="l"))
    })
    
    ## update filter choices upon change of data set
    shiny::observe({

        shiny::req(pgx$Y, pgx$samples)

        ## levels for sample filter
        levels = getLevels(pgx$Y)
        shiny::updateSelectInput(session, "data_samplefilter", choices=levels)

        grps <- pgx.getCategoricalPhenotypes(pgx$samples, min.ncat=2, max.ncat=999)
        grps <- sort(grps)
        selgrp <- grps[1]
        grps <- c("<ungrouped>",grps)
        if("group" %in% grps) selgrp = "group"
        if(nrow(pgx$samples)<=20) selgrp = "<ungrouped>"
        shiny::updateSelectInput(session,'data_groupby', choices=grps, selected=selgrp)
    })


      shiny::observeEvent({
          input$data_type
          pgx$X
          pgx$counts
      }, {

        if(input$data_type %in% c("counts","CPM")) {
            pp <- rownames(pgx$counts)
        } else {
            ## log2CPM
            pp <- rownames(pgx$X)
        }

        ## gene filter. 
        genes <- sort(pgx$genes[pp,]$gene_name)
        fc2 = rowMeans(pgx.getMetaFoldChangeMatrix(pgx)$fc**2)
        genes = intersect(names(sort(-fc2)),genes) ## most var gene??
        selgene <- genes[1]
        genes1 <- unique(c(selgene,sort(genes)))
        if(length(genes1)>1000) {
            genes1 <- c(sort(genes1[1:1000]),"(type SYMBOL for more genes...)",genes1[1001:length(genes1)])
        }
        shiny::updateSelectizeInput(session,'search_gene', choices=genes1, selected=selgene,
                                    ##options = list(maxOptions = 9999999),
                                    options = list(maxOptions = 1001),                                    
                                    server = TRUE)
    })

    last_search_gene <- reactiveVal()
    
    input_search_gene <- reactive({
        if( input$search_gene %in% c("(type SYMBOL for more genes...)","")) {
            gene1 <- last_search_gene() 
            return(gene1)
        }
        last_search_gene(input$search_gene)
        return(input$search_gene)
    })

    
    ##================================================================================
    ##=========================== MODULES ============================================
    ##================================================================================

    ## get selected samples after sample filtering
    selected_samples <- reactive({
        samples <- colnames(pgx$X)
        if(!is.null(input$data_samplefilter)) {
            samples <- selectSamplesFromSelectedLevels(pgx$Y, input$data_samplefilter)
        }
        samples
    })
    
    ## dbg("[***dataview_server] names.input = ",names(input))    
    dataview_module_geneinfo_server(
        "geneinfo", 
        r.gene  = reactive(input$search_gene)
    )

    ## first tab
    dataview_plot_tsne_server(
        "tsneplot",
        pgx,
        r.gene         = reactive(input$search_gene),
        r.samples      = selected_samples,
        r.data_type    = reactive(input$data_type),            
        r.data_groupby = reactive(input$data_groupby)
    )        

    dataview_plot_averagerank_server(
        "averagerankplot",
        pgx,
        r.gene         = reactive(input$search_gene),
        r.samples      = selected_samples,
        r.data_type    = reactive(input$data_type)            
    )            

    dataview_plot_correlation_server(
        "correlationplot",
        pgx,
        r.gene    = reactive(input$search_gene),
        r.samples = selected_samples
    )

    dataview_plot_tissue_server(
        "tissueplot",
        pgx,
        r.gene         = reactive(input$search_gene),
        r.data_type    = reactive(input$data_type)            
    )            

    dataview_plot_expression_server(
        "expressionplot",
        pgx,
        r.gene         = reactive(input$search_gene),
        r.samples      = selected_samples,
        r.data_type    = reactive(input$data_type),            
        r.data_groupby = reactive(input$data_groupby)
    )

    ## second tab
    dataview_plot_totalcounts_server("counts_total", input, getCountsTable)
    dataview_plot_boxplot_server("counts_boxplot", input, getCountsTable)
    dataview_plot_histogram_server("counts_histplot", input, getCountsTable)
    dataview_plot_abundance_server("counts_abundance", input, getCountsTable)
    dataview_plot_averagecounts_server("counts_averagecounts", input, getCountsTable)

    ## fourth tab
    dataview_plot_phenoheatmap_server("phenoheatmap", pgx, input)
    dataview_plot_phenoassociation_server("phenoassociation", pgx, input)    
    
    ##================================================================================
    ##========================= FUNCTIONS ============================================
    ##================================================================================
    
    getCountsTable <- shiny::reactive({
        shiny::req(pgx$X,pgx$Y,pgx$samples)

        shiny::validate(shiny::need("counts" %in% names(pgx), "no 'counts' in object."))
        subtt=NULL

        samples = colnames(pgx$X)
        samples <- selectSamplesFromSelectedLevels(pgx$Y, input$data_samplefilter)
        nsamples = length(samples)
        if("counts" %in% names(pgx)) {
            counts = pgx$counts[,samples,drop=FALSE]
        } else {
            cat("WARNING:: no counts table. estimating from X\n")
            counts = pmax(2**pgx$X - 1,0)
            k = grep("lib.size",colnames(pgx$samples))[1]
            if(length(k)>0) {
                libsize = pgx$samples[colnames(counts),k]
                libsize
                counts = t(t(counts) * libsize)
            }
            ##counts <- round(counts)
        }
        if(sum(is.na(counts))>0) {
            cat("WARNING:: plot counts: counts has missing values!\n")
        }

        ##if(input$data_sampling=="grouped") {
        grpvar <- input$data_groupby
        gr = pgx$Y[samples,grpvar]
        grps = sort(unique(gr))
        ##if(input$data_grouped && length(grps)>1 ) {
        if(input$data_groupby != "<ungrouped>" && length(grps)>1) {
            newx = c()
            for(g in grps) {
                mx = rowMeans(counts[,which(gr==g),drop=FALSE], na.rm=TRUE)
                ## mx = rowSums(counts[,which(gr==g),drop=FALSE], na.rm=TRUE)  ## SUM or MEAN???
                newx = cbind(newx, mx)
            }
            if(NCOL(newx)==1) newx <- matrix(newx,ncol=1)
            rownames(newx) = rownames(counts)
            colnames(newx) = grps
            counts = newx
        }

        ## if too many samples (like scRNA-seq do subsampling...)
        if(ncol(counts) > 500) {
            kk <- sample(ncol(counts),400,replace=TRUE)
            counts <- counts[,kk,drop=FALSE]
            subtt=c(subtt,"random subset")
        }
        colnames(counts) <- substring(colnames(counts),1,24)

        gset <- list()
        gg = pgx$genes[rownames(counts),]$gene_name
        tt = pgx$genes[rownames(counts),]$gene_title
        g1 <- gg[grep("^rpl|^rps",gg,ignore.case=TRUE)]
        g2 <- gg[grep("^mrpl|^mrps",gg,ignore.case=TRUE)]
        g3 <- gg[grep("^MT-",gg,ignore.case=TRUE)]
        g4 <- gg[grep("mitochondr",tt,ignore.case=TRUE)]
        gset[["Ribosomal (RPL/RPS)"]] = g1
        gset[["Mitochondrial ribosomal (MRPL/MRPS)"]] = g2
        gset[["Mitochondrial (MT)"]] = g3
        gset[["Other mitochondrial"]] = setdiff(g4,g3)
        jj <- grep("mitochondr|ribosom",names(FAMILIES),invert=TRUE,ignore.case=TRUE)
        gset.other <- lapply( FAMILIES[jj], function(x) setdiff(x, c(g1,g2,g3,g4)))
        gset <- c(gset, gset.other)
        gset <- gset[grep("<all>",names(gset),invert=TRUE)]
        gset <- gset[sapply(gset,length) > 10]

        ## Counts per samples, by category
        total.counts = Matrix::colSums(counts,na.rm=TRUE)
        summed.counts = t(sapply(gset, function(f)
            Matrix::colSums(counts[which(gg %in% f),,drop=FALSE], na.rm=TRUE)))
        avg.counts   = t(sapply(gset, function(f)
            Matrix::colMeans(counts[which(gg %in% f),,drop=FALSE], na.rm=TRUE)))
        prop.counts = 100 * t(t(summed.counts) / total.counts)

        head(sort(rowSums(prop.counts,na.rm=TRUE),decreasing=TRUE),10)
        head(sort(rowSums(avg.counts,na.rm=TRUE),decreasing=TRUE),10)
        ##jj <- order(-apply(avg.counts,1,sd,na.rm=TRUE))

        jj <- head(order(-rowSums(prop.counts,na.rm=TRUE)),6)
        prop.counts <- prop.counts[jj,,drop=FALSE]
        jj <- head(order(-rowSums(avg.counts,na.rm=TRUE)),6)
        avg.counts  <- avg.counts[jj,,drop=FALSE]
        sorting="no"
        if(sorting=="decr") {
            total.counts <- sort(total.counts, decreasing=TRUE)
            prop.counts  <- prop.counts[,order(-colMeans(prop.counts))]
            avg.counts   <- avg.counts[,order(-colMeans(avg.counts))]
            counts <- counts[,order(-colMeans(counts))]
        }
        if(sorting=="inc") {
            total.counts <- sort(total.counts, decreasing=FALSE)
            prop.counts  <- prop.counts[,order(colMeans(prop.counts))]
            avg.counts   <- avg.counts[,order(colMeans(avg.counts))]
            counts <- counts[,order(colMeans(counts))]
        }

        ss <- names(total.counts)
        prop.counts <- prop.counts[,ss,drop=FALSE]
        avg.counts <- avg.counts[,ss,drop=FALSE]
        counts <- counts[,ss,drop=FALSE]

        log2counts <- log2(1 + counts)
        ##log2counts[which(log2counts==0)] <- NA
        jj <- sample(nrow(counts),100,replace=TRUE)
        jj <- sample(nrow(counts),1000,replace=TRUE)

        ## create the plots
        par(mar=c(3,3,3,3), mgp=c(2.4,0.7,0), oma=c(1,1,1,1)*0.2 )
        cx1=1
        if(length(total.counts)>12) cx1=0.9
        if(length(total.counts)>30) cx1=0.8
        if(length(total.counts)>50) cx1=0.7
        if(length(total.counts)>80) cx1=0.6

        if(1) {
            names(total.counts) <- substring(names(total.counts),1,30)
            colnames(log2counts) <- substring(colnames(log2counts),1,30)
            colnames(prop.counts) <- substring(colnames(prop.counts),1,30)
            colnames(avg.counts) <- substring(colnames(avg.counts),1,30)
        }

        res = list(
            total.counts = total.counts,
            subtt = subtt,
            cx1 = cx1,
            log2counts = log2counts,
            jj = jj,
            prop.counts = prop.counts,
            avg.counts = avg.counts
        )
    })

    ##================================================================================
    ##===============================  TABLES ========================================
    ##================================================================================

    dropdown_search_gene='<code>Search gene</code>'
    menu_grouped='<code>grouped</code>'
    menu_options='<code>Options</code>'

    data_rawdataTable_text = paste0('Under the <strong>gene table </strong>, the average expression values of genes across the groups can be read. The samples (or cells) can be ungrouped by unclicking the ',menu_grouped, ' in the main <i>Options</i> to see the exact expression values per sample (or cell).', 'The genes in the table are ordered by the correlation (<b>rho</b> column) with respect to the gene selected by users from the ',dropdown_search_gene, ' setting. <b>SD</b> column reports the standard deviation of expression across samples (or cells).')

    data_rawdataTable.RENDER <- shiny::reactive({
        ## get current view of raw_counts

        shiny::req(pgx$X,pgx$Y,pgx$genes,pgx$model.parameters)
        shiny::req(input$data_type, input$data_groupby)

        dbg("[data_rawdataTable.RENDER] reacted")

        pp <- rownames(pgx$X)
        if(input$data_type=="counts") {
            ##x <- pgx$counts[pp,]
            x <- pgx$counts
        } else if(input$data_type=="CPM") {
            ##x <- 2**pgx$X
            ##x <- edgeR::cpm(pgx$counts[pp,], log=FALSE)
            x <- edgeR::cpm(pgx$counts, log=FALSE)
        } else {
            ## log2CPM
            x <- pgx$X
        }
        x0=x

        ##------------------ select samples
        dbg("[data_rawdataTable.RENDER] select samples")
        samples <- colnames(pgx$X)
        samples <- selectSamplesFromSelectedLevels(pgx$Y, input$data_samplefilter)
        samples <- intersect(colnames(x),samples)
        x <- x[,samples,drop=FALSE]

        gene = "CD3E"
        gene = "CCR6"
        gene = input_search_gene()
        if(is.null(gene) || gene=="" || is.na(gene)) return(NULL)
        ##xgene = sub(".*:","",rownames(x))

        ## Quickly (?) calculated correlation to selected gene
        dbg("[data_rawdataTable.RENDER] calculate rho")
        rho = sdx = avg = NULL
        if(input$data_type == "logCPM") {
            k=1
            logx <- pgx$X[rownames(x),]
            xgenes <- pgx$genes[rownames(x),"gene_name"]
            k <- which(xgenes==gene)
            rho = cor( t(logx[,samples]), logx[k,samples], use="pairwise")[,1]
            rho = round(rho[rownames(x)], digits=3)
            sdx = round(apply(logx[,samples],1,sd),digits=3)
        }
        avg <- round(rowMeans(x),digits=3)

        ##if(input$data_sampling=="grouped") {
        ##do.grouped <- input$data_grouped
        grpvar = "group"
        group <- NULL
        grpvar <- input$data_groupby

        if(grpvar %in% colnames(pgx$Y)) {
            group = pgx$Y[colnames(x),grpvar]
        }
        if(length(samples)>500 && grpvar=="<ungrouped>") {
            ##grpvar="group"
            group <- pgx$model.parameters$group
        }
        do.grouped <- (grpvar!="<ungrouped>")
        if(do.grouped && !is.null(group) ) {
            ##group = pgx$Y[colnames(x),grpvar]
            allgroups = sort(unique(group))
            newx = c()
            for(gr in allgroups) {
                mx = rowMeans(x[,which(group==gr),drop=FALSE],na.rm=TRUE)
                newx = cbind(newx, mx)
            }
            rownames(newx) = rownames(x)
            ##colnames(newx) = paste0("[",allgroups,"]")
            colnames(newx) = paste0("avg.",allgroups,"")
            x = newx
        }

        x = round( as.matrix(x), digits=3)
        x95 = quantile(as.vector(x0[which(x0>0)]),probs=0.95)
        x99 = quantile(as.vector(x0[which(x0>0)]),probs=0.99)

        if(NCOL(x)==0 || nrow(x)==0) return(NULL)

        dbg("[data_rawdataTable.RENDER] create dataframe")
        ##rownames(x) = sub(".*:","",rownames(x))
        xgenes <- pgx$genes[rownames(x),"gene_name"]
        gene.title <- GENE.TITLE[toupper(xgenes)]
        gene.title <- substring(gene.title,1,50)
        if(is.null(rho)) {
            x = data.frame( gene=xgenes, title=gene.title,
                           AVG=avg,
                           as.matrix(x), check.names=FALSE)
        } else {
            x = data.frame( gene=xgenes, title=gene.title,
                           rho=rho, SD=sdx, AVG=avg,
                           as.matrix(x), check.names=FALSE)
        }
        x = x[order(x$gene),,drop=FALSE]

        ## put selected gene always on top
        j1 <- which(x$gene==gene)
        jj <- c(j1, setdiff(1:nrow(x),j1))
        x <- x[jj,,drop=FALSE]

        if(ncol(x) > 100) {
            max.row <- 1e5 / ncol(x)
            max.row <- 100*ceiling(max.row/100)
            max.row
            x <- head(x, max.row)
        }
        numcols <- grep('gene|title',colnames(x),value=TRUE,invert=TRUE)

        DT::datatable( x, rownames=FALSE,
                      class = 'compact cell-border stripe hover',
                      extensions = c('Buttons','Scroller'),
                      selection = list(mode='single', target='row', selected=1),
                      options=list(
                          dom = 'lfrtip',
                          ##pageLength = 60,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scroller=TRUE, scrollX = TRUE, scrollY = tabH,
                          deferRender=TRUE
                      )  ## end of options.list
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
                DT::formatStyle(numcols,
                                background = DT::styleColorBar(c(0,x99), 'lightblue'),
                                ##background = color_from_middle(x99, 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')
    }) %>%
    bindCache(input$search_gene, input$data_type, input$data_groupby)


    data_rawdataTable <- shiny::callModule(
        tableModule, "data_rawdataTable",
        func = data_rawdataTable.RENDER,
        title = "Gene expression table",
        filename = "counts.csv",
        info.text = data_rawdataTable_text
    )

    ##================================================================================
    ##================================= Samples ======================================
    ##================================================================================

    data_sampleTable.RENDER <- shiny::reactive({
        ## get current view of raw_counts
        shiny::req(pgx$Y,pgx$samples)
        
        dt <- NULL
        samples <- selectSamplesFromSelectedLevels(pgx$Y, input$data_samplefilter)
        dt <- pgx$samples[samples,,drop=FALSE]
        DT::datatable( dt,
                      class = 'compact cell-border stripe hover',
                      rownames = TRUE,
                      extensions = c('Buttons','Scroller'),
                      selection = list(mode='single', target='row', selected=1),
                      options=list(
                          dom = 'lfrtip',
                          scroller=TRUE, scrollX = TRUE, scrollY = 190,
                          deferRender=TRUE
                      )) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')         
    }) %>%
    bindCache(input$data_samplefilter)
    
    data_sampleTable_info = "<b>Sample information table.</b> Phenotype information about the samples. Phenotype variables starting with a 'dot' (e.g. '.cell cycle' and '.gender' ) have been estimated from the data."

    data_sampleTable <- shiny::callModule(
        tableModule, "data_sampleTable", label="c",
        func = data_sampleTable.RENDER,
        func2 = data_sampleTable.RENDER,
        title = "Sample information",
        filename = "samples.csv",
        info.text = data_sampleTable_info,
        height = c(280,750), width=c('auto',1280)
    )

    ##================================================================================
    ##================================= CONTRASTS ====================================
    ##================================================================================

    data_contrastTable.RENDER <- shiny::reactive({
        ## get current view of raw_counts

        shiny::req(pgx$Y,pgx$model.parameters)

        dbg("[data_contrastTable.RENDER] reacted")

        ##if(is.null(input$data_samplefilter)) return(NULL)
        dt <- NULL
        samples <- selectSamplesFromSelectedLevels(pgx$Y, input$data_samplefilter)
        names(pgx$model.parameters)
        if(input$data_ctbygroup=="group") {
            ct <- pgx$model.parameters$contr.matrix
            ##kk <- which(rownames(ct) %in% pgx$samples[samples,"group"])
            kk <- which(rownames(ct) %in% pgx$model.parameters$group[samples])
            dt <- ct[kk,,drop=FALSE]
        } else {
            dt <- pgx$model.parameters$exp.matrix[samples,,drop=FALSE]
        }
        dt <- sign(dt)
        colnames(dt) <- sub("[_. ]vs[_. ]","\nvs ",colnames(dt))
        dt[dt==0] <- NA

        DT::datatable( dt,
                      class = 'compact cell-border stripe hover',
                      rownames = TRUE,
                      extensions = c('Buttons','Scroller'),
                      selection = list(mode='single', target='row', selected=1),
                      options=list(
                          dom = 'lfrtip',
                          scroller=TRUE, scrollX = TRUE, scrollY = tabH,
                          deferRender=TRUE,
                          autoWidth = TRUE
                      )) %>%
            DT::formatStyle(0, target='row', fontSize='12px', lineHeight='70%') %>%
                DT::formatStyle(colnames(dt),
                                background = color_from_middle(c(-1,1), 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')

    })

    data_contrastTable_info = "<b>Contrast table.</b> Table summarizing the contrasts of all comparisons. Here, you can check which samples belong to which groups for the different comparisons. Non-zero entries '+1' and '-1' correspond to the group of interest and control group, respectively. Zero or empty entries denote samples not use for that comparison."


    data_contrastTable_opts = shiny::tagList(
        withTooltip( shiny::radioButtons(ns('data_ctbygroup'),
                             "Show by:", choices=c("sample","group")),
               "Show contrasts by group or by samples.",
               placement="right", options = list(container = "body"))
    )

    data_contrastTable <- shiny::callModule(
        tableModule, "data_contrastTable",
        func = data_contrastTable.RENDER,
        options = data_contrastTable_opts,
        title = "Contrast table",
        filename = "contrasts.csv",
        info.text = data_contrastTable_info
    )

    ##================================================================================
    ## Resource info (dev)
    ##================================================================================

    datatable_timings.RENDER <- shiny::reactive({

        shiny::req(pgx$timings)

        dbg("[datatable_timings.RENDER] reacted")

        ##if(is.null(pgx$timings)) return(NULL)
        D <- data.frame()
        if(!is.null(pgx$timings)) {
            D <- round(pgx$timings[,1:3],digits=3)
            D <- apply(D, 2, function(x) tapply(x, rownames(D), sum))
            catg <- gsub("^\\[|\\].*","",rownames(D))
            metd <- gsub("^.*\\]","",rownames(D))
            D <- data.frame(category=catg, method=metd, D)
        }
        DT::datatable( D, rownames=FALSE,
                      options = list(dom='tp', pageLength = 100),
                      class = 'compact cell-border stripe hover' ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')
    })

    datatable_timings_text = 'The <b>timings</b> table reports more detailed information about the object dimensions, object sizes and execution times of the methods.'

    datatable_timings <- shiny::callModule(
        tableModule, "datatable_timings",
        func = datatable_timings.RENDER,
        info.text = datatable_timings_text,
        options = NULL, title='Timings'
    )

    datatable_objectdims.RENDER <- shiny::reactive({

        shiny::req(pgx$X)

        dims1 <- lapply( pgx, dim)
        lens <- sapply( pgx, length)
        dims2 <- t(sapply( pgx[which(!sapply(dims1,is.null)) ], dim))
        kk <- which(sapply(dims1,is.null))
        dims2 <- rbind(dims2, cbind(lens[kk],0))
        colnames(dims2) = c("nrows","ncols")
        D = data.frame( object=rownames(dims2), dims2, check.names=FALSE)
        DT::datatable( D, rownames=FALSE,
                      options = list(dom='t', pageLength = 50),
                      class = 'compact cell-border stripe hover') %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')
    })

    datatable_objectdims_text = 'This table provides details about the data dimensions of objects.'

    datatable_objectdims <- shiny::callModule(
        tableModule, "datatable_objectdims",
        func = datatable_objectdims.RENDER,
        info.text = datatable_objectdims_text,
        options = NULL, title='Object dimensions'
    )

    datatable_objectsize.RENDER <- shiny::reactive({

        shiny::req(pgx$name)

        objsize <- sapply(pgx,object.size)
        objsize <- round( objsize/1e6, digits=2)
        D = data.frame( object=names(pgx), "size.Mb"=objsize, check.names=FALSE)
        DT::datatable( D, rownames=FALSE,
                      options = list(dom='t', pageLength = 50),
                      class = 'compact cell-border stripe hover') %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')
    })

    datatable_objectsize_text = "This table provides information about  about the memory sizes of objects"

    datatable_objectsize <- shiny::callModule(
        tableModule, "datatable_objectsize",
        func = datatable_objectsize.RENDER,
        options = NULL, title='Object sizes',
        info.text = datatable_objectsize_text
    )
  })
}
