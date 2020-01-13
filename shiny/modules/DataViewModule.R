
DataViewInputs <- function(id) {
    ns <- NS(id)  ## namespace
    tagList(
        uiOutput(ns("description")),
        uiOutput(ns("inputsUI"))
    )
}

DataViewUI <- function(id) {
    ns <- NS(id)  ## namespace
    tabsetPanel(
        id = ns("tabs"),
        tabPanel("Plots",uiOutput(ns("plotsUI"))),
        tabPanel("Counts",uiOutput(ns("countsUI"))),
        tabPanel("Gene table",uiOutput(ns("genetableUI"))),
        tabPanel("Samples",uiOutput(ns("sampletableUI"))),
        tabPanel("Contrasts",uiOutput(ns("contrasttableUI"))),        
        tabPanel("Resource info",uiOutput(ns("resourceinfoUI")))
    )
}

DataViewModule <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE
    inputData <- env[["load"]][["inputData"]]
    
    rowH = 355  ## row height of panels
    imgH = 315  ## height of images
    fullH = 750 ## full height of panel
    tabH = 600  ## height of tables
    
    description = "<b>DataView.</b> Information and descriptive statistics to quickly lookup a gene, check the total counts, or view the data tables."
    output$description <- renderUI(HTML(description))
    
    ##----------------------------------------------------------------------
    ## More Info (pop up window)
    ##----------------------------------------------------------------------
    dropdown_search_gene='<code>Search gene</code>'
    menu_grouped='<code>Group by</code>'
    
    data_infotext =paste0(
        'The <strong>DataView module</strong> provides information and visualisations of the dataset to quickly lookup a gene, check the counts, or view the data tables. 

<br><br>The <strong>Plots</strong> panel displays figures related to the expression level of the selected gene, correlation, and average expression ranking within the dataset. More information about the gene and hyperlinks to external databases are provided. Furthermore, it displays the correlation and tissue expression for a selected gene in external reference datasets. In the <strong>Counts</strong> panel, the total number of counts (abundance) per sample and their distribution among the samples are displayed. This is most useful to check the technical quality of the dataset, such as total read counts or abundance of ribosomal genes. In <strong>Gene Table</strong> panel, the exact expression values across the samples can be looked up, where genes are ordered by the correlation with respect to the first gene. Gene-wise average expression of a phenotype sample grouping is also presented in this table. In the <strong>Samples</strong> panel, more complete information about samples can be found. Finally, the <strong>Contrasts</strong> panel, shows information about the phenotype comparisons.

<br><br><br><br>
<center><iframe width="500" height="333" src="https://www.youtube.com/embed/watch?v=IrqU-z9UDqc&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=1" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>

')

    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================

    require(htmltools)
    ## data set parameters
    datatypes <- c("CPM","logCPM")
    datatypes <- c("counts","CPM","logCPM")

    output$inputsUI <- renderUI({
        ui <- tagList(
            tipify( actionLink(ns("data_info"), "Tutorial", icon = icon("youtube")),
                   "Show more information about this module."),
            hr(), br(), 
            ##textInput("search_gene","Filter genes", value="")
            tipify( selectInput(ns("search_gene"),"Search gene:", choices=NULL),
                   "Select a gene of interest for the analysis.", placement="top"),
            tipify( selectInput(ns("data_samplefilter"),"Filter samples:",
                                choices=NULL, multiple=TRUE),
                   "Filter the relevant samples for the analysis.", placement="top"),
            tipify( selectInput(ns('data_groupby'),'Group by:', choices=NULL),
                   "Select phenotype for grouping the samples.", placement="top"),
            br(),
            tipify( actionLink(ns("data_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                   "Toggle advanced options.", placement="top"),
            br(),br(),
            conditionalPanel(
                "input.data_options % 2 == 1", ns=ns,
                tipify( radioButtons(ns('data_type'),'Data type:',
                                     choices=datatypes, selected="logCPM", inline=TRUE),
                       "Choose an input data type for the analysis.", placement="bottom")
            )
        )
        ui
    })
    outputOptions(output, "inputsUI", suspendWhenHidden=FALSE) ## important!!!
    
    ##================================================================================
    ##========================= OBSERVE ==============================================
    ##================================================================================

    ## ------- observe functions -----------
    observeEvent( input$data_info, {
        showModal(modalDialog(
            title = HTML("<strong>Data View Module</strong>"),
            HTML(data_infotext),
            easyClose = TRUE, size="l" ))
    })

    ## update filter choices upon change of data set 
    observe({
        ngs <- inputData()
        req(ngs)

        ## levels for sample filter
        levels = getLevels(ngs$Y)
        updateSelectInput(session, "data_samplefilter", choices=levels)
        genes <- sort(ngs$genes[rownames(ngs$X),]$gene_name)
        sel = genes[1]  ## most var gene
        updateSelectInput(session,'search_gene', choices=genes, selected=sel)
        
        grps <- pgx.getCategoricalPhenotypes(ngs$samples, min.ncat=2, max.ncat=20)
        grps <- c("<ungrouped>",sort(grps))
        selgrp <- grps[1]
        if("group" %in% grps) selgrp = "group"
        if(nrow(ngs$samples)<=20) selgrp = "<ungrouped>"
        updateSelectInput(session,'data_groupby', choices=grps, selected=selgrp)
        
    })
    

    ##================================================================================
    ##========================= FUNCTIONS ============================================
    ##================================================================================


    ##----------------------------------------------------------------------
    ##                     Info messsages
    ##----------------------------------------------------------------------
    
    data_genePlots_tsne_text=paste0('<b>T-SNE clustering</b> of samples (or cells) colored by an expression of the gene selected in the ',dropdown_search_gene, ' dropdown menu. The red color represents an over-expression of the selected gene across samples (or cells).')
    data_genePlots_barplot_text=paste0('Expression barplot of grouped samples (or cells) for the gene selected in the ',dropdown_search_gene,'. Samples (or cells) in the barplot can be ungrouped by setting the ',menu_grouped, ' under the main <i>Options</i>.')
    data_genePlots_correlationplot_text=paste0('Barplot of the top positively and negatively correlated genes with the selected gene. Absolute expression levels of genes are colored in the barplot, where the low and high expressions range between the light and dark colors, respectively.')
    data_genePlots_averageRankPlot_text=paste0('Ranking of the average expression of the selected gene.')
    data_geneInfo_text = paste0('To find out more information from the literature, hyperlinks are provide to connect the selected gene to public databases, including ', a_OMIM,', ', a_KEGG, ' and ',a_GO,'.')
    data_corplot_text = paste0('Top cumulative positively and negatively correlated genes with the selected gene in the current dataset as well as in public datasets such as ',a_ImmProt,' and ',a_HPA,'. The correlations of genes are colored by dataset.')
    data_tissueplot_text = paste0('Tissue expression for the selected gene in the tissue expression ',a_GTEx,' dataset. Colors corresponds to "tissue clusters" as computed by unsupervised clustering.')
    

    ##----------------------------------------------------------------------
    ##                     Average Rank plot
    ##----------------------------------------------------------------------
    MARGINS1 = c(7,3.5,2,1)

    data_genePlots_averageRankPlot.RENDER %<a-% reactive({
        require(RColorBrewer)
        
        ngs <- inputData()
        req(ngs)
        
        dbg("[data_genePlots_averageRankPlot.RENDER] reacted")

        gene = "KCNN4"
        gene = ngs$genes$gene_name[1]
        if(!is.null(input$search_gene) && input$search_gene!="") gene <- input$search_gene
        samples = colnames(ngs$X)
        if(!is.null(input$data_samplefilter)) {
            samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        }
        nsamples = length(samples)
        
        if(input$data_type=="counts") {
            mean.fc <- sort(rowMeans(ngs$counts[,samples,drop=FALSE]),decreasing=TRUE)
            ylab="expression (counts)"
        }
        if(input$data_type=="CPM") {
            ##cpm <- edgeR::cpm(ngs$counts[,samples])
            cpm = 2**ngs$X[,samples,drop=FALSE]
            mean.fc <- sort(rowMeans(cpm),decreasing=TRUE)
            ylab="expression (CPM)"
        }
        if(input$data_type=="logCPM") {
            mean.fc <- sort(rowMeans(ngs$X[,samples,drop=FALSE]),decreasing=TRUE)
            ylab="expression (log2CPM)"
        }
        
        j <- which(sub(".*:","",names(mean.fc))==gene)
        ##j <- which(ngs$genes$gene_name==gene)
        
        mar = MARGINS1
        par(mar=mar, mgp=c(2.1,0.8,0))
        ##MARGINS1
        plot( mean.fc, type="h", lwd=0.4,
             ## col="#3380CC11", 
             col="#bbd4ee", 
                                        #main="average rank", cex.main=1.2,
             ylab=ylab, xlab="ordered genes", xaxt="n")
        points( j, mean.fc[j], type="h", lwd=2, col="black")
        text( j, mean.fc[j], gene, pos=3, cex=0.9)

        dbg("[data_genePlots_averageRankPlot.RENDER] done")

    })

    ##data_genePlots_averageRankPlot_module <- plotModule(
    ##id="data_genePlots_averageRankPlot", ns=ns,
    callModule(
        plotModule, id="data_genePlots_averageRankPlot",
        func = data_genePlots_averageRankPlot.RENDER,
        func2 = data_genePlots_averageRankPlot.RENDER,
        info.text = data_genePlots_averageRankPlot_text,
        height = imgH, ## width = '100%',
        pdf.width=6, pdf.height=6, 
        label="d", title="Average rank"
    )
    ##output <- attachModule(output, data_genePlots_averageRankPlot_module) 

    ##----------------------------------------------------------------------
    ##                     Correlation plot
    ##----------------------------------------------------------------------
    
    data_genePlots_correlationplot_data <- reactive({
        require(RColorBrewer)
        
        ngs <- inputData()
        req(ngs)

        dbg("[data_genePlots_correlationplot_data()] reacted")
        
        gene = "KCNN4"
        gene = ngs$genes$gene_name[1]
        if(!is.null(input$search_gene) && input$search_gene!="") gene <- input$search_gene
        samples = colnames(ngs$X)
        if(!is.null(input$data_samplefilter)) {
            samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        }
        nsamples = length(samples)
        
        grp = factor(ngs$Y[samples,"group"])
        klr0 = rep(brewer.pal(8,"Set2"),99)
        klr0 = COLORS
        klr = klr0[factor(ngs$Y[samples,"group"])]
        
        ## precompute
        pp=rownames(ngs$genes)[1]
        pp <- rownames(ngs$genes)[match(gene,ngs$genes$gene_name)]
        
        gx = NULL
        ylab = NULL
        if(input$data_type=="counts") {
            gx = ngs$counts[pp,samples]
            ylab="expression (counts)"
        } else if(input$data_type=="CPM") {
            gx = 2**ngs$X[pp,samples]
            ylab="expression (CPM)"
        } else if(input$data_type=="logCPM") {
            gx = ngs$X[pp,samples]
            ylab="expression (log2CPM)"
        }
        
        mar = MARGINS1
        par(mar=mar, mgp=c(2.1,0.8,0))
        
        ## corr always in log.scale and restricted to selected samples subset
        ## should match exactly the rawtable!!
        rho = cor(t(ngs$X[,samples]), ngs$X[pp,samples], use="pairwise")[,1]
        rho[is.na(rho)] <- 0
        jj = head(order(-abs(rho)),30)
        jj <- c( head(order(rho),15), head(order(-rho),15))
        jj <- jj[order(-rho[jj])]
        top.rho = rho[jj]
        
        gx1 <- sqrt(rowSums(ngs$X[names(top.rho),samples]**2,na.rm=TRUE))
        gx1 <- (gx1 / max(gx1))
        ## klr1 <- rev(grey.colors(16,start=0.2,end=0.9,gamma=0.25))[1+round(15*gx1) ]
        ## klr1[which(is.na(klr1))] <- "#DDDDDD"
        ##------ Color test ------#
        klr1 <- rev(colorRampPalette(c(rgb(0.2,0.5,0.8,0.8), rgb(0.2,0.5,0.8,0.1)), alpha = TRUE)(16))[1+round(15*gx1) ]
        klr1[which(is.na(klr1))] <- rgb(0.2,0.5,0.8,0.1)
        ##------------------------#
        names(top.rho) = sub(".*:","",names(top.rho))
        offset = min(top.rho)*0.95
        offset = 0

        dbg("[data_genePlots_correlationplot_data()] done")
        
        res = list(top.rho=top.rho, offset=offset, klr1=klr1)
        return(res)
    })

    data_genePlots_correlationplot.RENDER %<a-% reactive({

        dbg("[data_genePlots_correlationplot.RENDER] reacted")

        res = data_genePlots_correlationplot_data()
        ngs <- inputData()
        if(is.null(ngs) | is.null(res)) return(NULL)

        gene = ngs$genes$gene_name[1]
        if(!is.null(input$search_gene) && input$search_gene!="") gene <- input$search_gene
        barplot(res$top.rho - res$offset, col=res$klr1, ## horiz=TRUE,
                las=3,#main=paste("top correlated genes\nwith",gene),
                main=paste(gene),
                offset = res$offset, ylab="correlation (r)",
                ## names.arg=rep(NA,length(top.rho)),
                cex.names=0.65, cex.main=1, col.main="#7f7f7f", border=NA)
        ##text( (1:length(top.rho) - 0.5)*1.2, offset, names(top.rho),
        ##col="black", cex=0.75, srt=90, pos=3, offset=0.4, font=1)
        legend("topright", legend=c("expr high","expr low"),
                                        #fill=c("grey30","grey80"),
               fill=c("#3380CCCC","#3380CC40"),
               cex=0.85, y.intersp=0.85)

        dbg("[data_genePlots_correlationplot.RENDER] done")
        
    })
    
    ##data_genePlots_correlationplot_module <- plotModule(
    ##    id="data_genePlots_correlationplot", ns=ns,
    callModule(
        plotModule, "data_genePlots_correlationplot",
        func = data_genePlots_correlationplot.RENDER,
        func2 = data_genePlots_correlationplot.RENDER,
        info.text=data_genePlots_correlationplot_text,
        height = imgH, pdf.width=6, pdf.height=6,
        label="c", title="Top correlated genes"
    )
    ##output <- attachModule(output, data_genePlots_correlationplot_module) 

    ##----------------------------------------------------------------------
    ##                     Bar/box plot
    ##---------------------------------------------------------------------- 
    
    data_genePlots_barplot.RENDER %<a-% reactive({
        require(RColorBrewer)
        
        cat("[dataview] data_genePlots_barplot.RENDER reacted\n")
        
        ngs <- inputData()
        req(ngs)
        req(input$data_groupby,input$search_gene,input$data_type)
                
        gene = "KCNN4"
        gene = ngs$genes$gene_name[1]
        if(!is.null(input$search_gene) && input$search_gene!="") gene <- input$search_gene
        samples = colnames(ngs$X)
        if(!is.null(input$data_samplefilter)) {
            samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        }
        nsamples = length(samples)
        
        grpvar <- input$data_groupby
        ##grp  = factor(ngs$Y[samples,grpvar])
        grp  = factor(as.character(ngs$Y[samples,grpvar]))
        klr0 = COLORS
        klr  = klr0[as.integer(grp)]
        
        ## precompute
        pp=rownames(ngs$genes)[1]
        pp <- rownames(ngs$genes)[match(gene,ngs$genes$gene_name)]
        
        gx = NULL
        ylab = NULL
        if(input$data_type=="counts") {
            gx = ngs$counts[pp,samples]
            ylab="expression (counts)"
        } else if(input$data_type=="CPM") {
            gx = 2**ngs$X[pp,samples]
            ylab="expression (CPM)"
        } else if(input$data_type=="logCPM") {
            gx = ngs$X[pp,samples]
            ylab="expression (log2CPM)"
        }
                
        mar=MARGINS1
        par(mar=mar, mgp=c(2.1,0.8,0))

        BLUE = col=rgb(0.2,0.5,0.8,0.8)
        bee.cex = ifelse(length(gx)>500,0.1,0.2)
        bee.cex = c(0.3,0.1,0.05)[cut(length(gx),c(0,100,500,99999))]

        ##if(input$data_grouped) {
        if(input$data_groupby != "<ungrouped>") {
            ngrp <- length(unique(grp))
            cx1 = ifelse( ngrp < 10, 1, 0.8)
            cx1 = ifelse( ngrp > 20, 0.6, cx1)
                                        #cx1 = ifelse( ngrp > 10, 0.6, 0.9)
            gx.b3plot( gx, grp, las=3, main=gene, ylab=ylab, 
                      cex.main=1, col.main="#7f7f7f",
                      bar=TRUE, border=NA, ## bee = ifelse(length(gx) < 500,TRUE,FALSE), 
                      bee.cex=bee.cex, ## sig.stars=TRUE, max.stars=5,
                      xlab="", names.cex=cx1, srt=45,
                      ## col=klr0[ii], 
                      col=rgb(0.2,0.5,0.8,0.8)) #col=rgb(0.2,0.5,0.9,0.6)) 
        }  else {
            jj <- 1:length(gx)
            
            sorting="no"
            if(sorting == "decr")  jj <- order(-gx)
            if(sorting == "inc")  jj <- order(gx)
            tt=""
            barplot(gx[jj], col=BLUE, ##col=klr[jj],
                    las=3, cex.names=0.85,
                    ylab=ylab, xlab=tt,
                    main=gene, cex.main=1, col.main="#7f7f7f", border=NA,
                    names.arg=rep(NA,length(gx)) )
            if(length(gx)<100) {
                cx1 = ifelse(length(gx) > 20, 0.8, 0.9)
                cx1 = ifelse(length(gx) > 40, 0.6, cx1)
                cx1 = ifelse(length(gx) < 10, 1, cx1)
                text((1:length(gx)-0.5)*1.2, -0.04*max(gx), names(gx)[jj], 
                     las=3, cex=cx1, pos=2, adj=0, offset=0, srt=60, xpd=TRUE)
            }
        }
    })

    ##data_genePlots_barplot_module <- plotModule(
    ##    id="data_genePlots_barplot", ns=ns,
    callModule(
        plotModule, "data_genePlots_barplot",    
        func = data_genePlots_barplot.RENDER,
        func2 = data_genePlots_barplot.RENDER,
        info.text=data_genePlots_barplot_text,
        height=imgH, pdf.width=6, pdf.height=6, ## height="400px",    
        label="b", title="Abundance/expression barplot"
    )
    ##output <- attachModule(output, data_genePlots_barplot_module)

    ##----------------------------------------------------------------------
    ## t-SNE
    ##----------------------------------------------------------------------
    
    data_genePlots_tsne.RENDER %<a-% reactive({
        require(RColorBrewer)
        
        ngs <- inputData()
        req(ngs)

        dbg("[data_genePlots_tsne.RENDER] reacted")
        
        gene = "KCNN4"
        gene = ngs$genes$gene_name[1]
        if(!is.null(input$search_gene) && input$search_gene!="") gene <- input$search_gene
        samples = colnames(ngs$X)
        if(!is.null(input$data_samplefilter)) {
            samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        }
        nsamples = length(samples)
        
        grp = factor(ngs$Y[samples,"group"])    
        klr0 = rep(brewer.pal(8,"Set2"),99)
        klr0 = COLORS
        klr = klr0[factor(ngs$Y[samples,"group"])]
        
        ## precompute
        pp=rownames(ngs$genes)[1]
        pp <- rownames(ngs$genes)[match(gene,ngs$genes$gene_name)]
        
        gx = NULL
        ylab = NULL
        if(input$data_type=="counts") {
            gx = ngs$counts[pp,samples]
            ylab="expression (counts)"
        } else if(input$data_type=="CPM") {
            gx = 2**ngs$X[pp,samples]
            ylab="expression (CPM)"
        } else if(input$data_type=="logCPM") {
            gx = ngs$X[pp,samples]
            ylab="expression (log2CPM)"
        }
        
        ##par(mar=c(12,2,2,1), mgp=c(2.1,0.8,0), oma=c(3,0.5,1.5,0.3))
        require(RColorBrewer)
        pos <- ngs$tsne2d[samples,]
        
        cex1 <- 1.8*c(1.6,1.0,0.6,0.3)[cut(nrow(pos),breaks=c(-1,40,200,1000,1e10))]    
        klrpal = colorRampPalette(c("blue3", "aliceblue", "grey85", "lavenderblush", "red3"))(16)
        klrpal = colorRampPalette(c("grey80", "grey50", "red3"))(16)
        
        fc1 <- tanh(0.99 * scale(gx)[,1])
        fc1 <- tanh(0.99 * scale(gx,center=FALSE)[,1])
        ##fc1 <- tanh(0.99 * gx/sd(gx))        
        fc2 <- (fc1 - min(fc1))
        klr1 = klrpal[1 + round(15*fc2/max(abs(fc2)))]
        klr1 = paste0(col2hex(klr1),"88")
        
        ##par(mar=c(8,2,2.2,1), mgp=c(1,0.5,0))
        par(mar=c(2.3,2.3,2,2), mgp=c(0.9,0.1,0))
        jj2 <- order(abs(fc1))
        plot( pos[jj2,], pch=20, cex=cex1, col=klr1[jj2], fg = gray(0.6), bty = "o",
             xaxt='n', yaxt='n', xlab="tSNE1", ylab="tSNE2")
        
        ## determine how to do grouping for group labels
        grp <- factor(ngs$samples[samples,]$group)
        ngrp <- length(unique(grp))
        
        if("cell.type" %in% colnames(ngs$samples) && ngrp>20) {
            grp <- ngs$samples[samples,]$cell.type
            dbg("[DataView] change to cell.type grp : head.grp=",head(grp))
        }
        
        cex2 = ifelse(nrow(pos) < 50, 1.5, 1.1)
        cex2 = ifelse(nrow(pos) > 200, 0.8, cex2)
        if(0 && input$pr_labelmode=="legend") {
            legend("bottomright", legend=levels(grp), fill=klrpal,
                   cex=cex2, y.intersp=0.8, bg="white")
        } else {
            ##grp.pos <- apply(pos,2,function(x) tapply(x,grp,mean))
            grp.pos <- apply(pos,2,function(x) tapply(x,grp,median))
            if(length(unique(grp))==1) {
                grp.pos <- matrix(grp.pos,ncol=2)
                rownames(grp.pos) <- unique(grp)
            }
            labels = rownames(grp.pos)
            cex3 <- c(1.4,1.2,1,0.8)[cut(length(labels),breaks=c(-1,5,10,20,999))]
            boxes = sapply(nchar(labels),function(n) paste(rep("\u2588",n),collapse=""))
            text( grp.pos, labels=boxes, cex=0.9*cex3, col="#CCCCCC88")
            text( grp.pos, labels=labels, font=2, cex=0.9*cex3, col="black")
            ##text( grp.pos[,], labels=rownames(grp.pos), font=2, cex=cex1**0.5)
        }

        dbg("[data_genePlots_tsne.RENDER] done")
        
    })

    ##data_genePlots_tsne_module <- plotModule(
    ##    id="data_genePlots_tsne", ns=ns,
    callModule(
        plotModule, "data_genePlots_tsne",
        func = data_genePlots_tsne.RENDER,
        func2 = data_genePlots_tsne.RENDER,
        info.text = data_genePlots_tsne_text,
        height = imgH, pdf.width = 6, pdf.height = 6,
        label = "a", title= "t-SNE clustering"
    )
    ##output <- attachModule(output, data_genePlots_tsne_module)
    
    ##----------------------------------------------------------------------
    ##  Tissue expression plot
    ##----------------------------------------------------------------------
    
    data_tissueplot.RENDER  %<a-% reactive({
        
        ngs <- inputData()
        req(ngs)
        if(is.null(input$data_type)) return(NULL)

        dbg("[data_tissueplot.RENDER] reacted")
        
        gene <- input$search_gene
        pp <- rownames(ngs$genes)[match(gene,ngs$genes$gene_name)]    
        hgnc.gene = toupper(as.character(ngs$genes[pp,"gene_name"]))
        
        require(RColorBrewer)
        par(mar=c(6,4,1,1), mgp=c(2.2,0.8,0))
        mar=MARGINS1
        par(mar=mar, mgp=c(1.5,0.5,0))
        
        if( hgnc.gene %in% rownames(TISSUE)) {
            tx = TISSUE[hgnc.gene,]
            grp = TISSUE.grp[names(tx)]
            tissue.klr = COLORS[grp]
            ylab="expression (TPM)"
            if(input$data_type=="logCPM") {
                ylab = "expression (log2TPM)"
                tx = log(1 + tx)
            }
            jj <- 1:length(tx)
            sorting="no"
            if(sorting=="decr") jj <- order(-tx)
            if(sorting=="inc") jj <- order(tx)
            
            barplot(tx[jj], las=3, main=hgnc.gene, cex.main=1, col.main="#7f7f7f",
                    col = tissue.klr[jj], border=NA,ylab=ylab, cex.names=0.9,
                    names.arg=rep(NA,length(tx)))
                                        #title(main="tissue expression", line=-0.3, cex.main=1.1)
            text((1:length(tx)-0.5)*1.2, -0.04*max(tx), names(tx)[jj], las=3,
                 cex=0.85, pos=2, adj=0, offset=0, srt=55, xpd=TRUE)
            
        } else {
            frame()
        }

        dbg("[data_tissueplot.RENDER] done")
        
    })
    
    ##data_tissueplot_module <- plotModule(
    ##    id="data_tissueplot", ns=ns,
    callModule(
        plotModule, "data_tissueplot",    
        func = data_tissueplot.RENDER,
        func2 = data_tissueplot.RENDER,
        info.text = data_tissueplot_text,
        height = imgH, pdf.width=9, pdf.height=6,
        label="g", title="Tissue expression"
    )
    ##output <- attachModule(output, data_tissueplot_module) 
    
    ##----------------------------------------------------------------------
    ##  Cumulative correlation plot (stacked)
    ##----------------------------------------------------------------------
    
    data_corplot_data <- reactive({
        require(RColorBrewer)
        ngs <- inputData()
        req(ngs)	
        if(is.null(input$data_type)) return(NULL)

        dbg("[data_corplot_data()] reacted")
        
        samples=colnames(ngs$X);gene="CD4"
        samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        gene <- input$search_gene
        if(is.null(gene)) return(NULL)    
        
        ## corr always in log.scale and restricted to selected samples subset
        zx <- ngs$X
        grp <- ngs$samples$group
        dim(zx)
        if( FALSE && length(grp) >= 5 && ncol(zx) > 50) {
            ## TOO SLOW!!! should do pre-computed???
            zx <- t( apply(ngs$X, 1, function(x) tapply(x,grp,mean)))
        }
        head(zx)
        rownames(zx) <- toupper(sub(".*:","",rownames(zx)))  ## NEED RETHINK!
        
        xref <- list("this_data"=zx, HPA_tissue=as.matrix(TISSUE),
                     ImmProt=as.matrix(IMMPROT))
        gene0 <- toupper(gene)  ## uppercase mouse
        R <- pgx.getGeneCorrelation(gene0, xref=xref)    
        if(is.null(R)) return(NULL)
        
        rho.genes = as.character(ngs$genes$gene_name)
        if("hgnc_symbol" %in% colnames(ngs$genes)) {
            rho.genes = as.character(ngs$genes$hgnc_symbol)
        }
        R <- R[match(rho.genes,rownames(R)),,drop=FALSE]
        rownames(R) <- rho.genes
        
        ## get top correlated genes
        ##jj = head(order(rowSums(R),decreasing=FALSE),35)
        rsum <- rowSums(R,na.rm=TRUE)
        jj = head(order(abs(rsum),decreasing=TRUE),35)
        jj = head(order(-abs(rsum)),30)
        jj <- c( head(order(rsum),20), head(order(-rsum),20))
        jj <- jj[order(-rsum[jj])]
        head(rsum[jj])
        Rtop = R[jj,,drop=FALSE]
        rownames(Rtop) = sub(".*:","",rownames(Rtop))
        offset = min(Rtop, na.rm=TRUE)*0.95
        offset=0
        klr <- grey.colors(ncol(Rtop),start=0.3,end=0.7)
        
        ## --- color test -----##
        klr <- grey.colors(ncol(Rtop),start=0.3,end=0.7)
        klr <- colorRampPalette(c(rgb(0.2,0.5,0.8,0.8), rgb(0.2,0.5,0.8,0.2)), alpha = TRUE)(ncol(Rtop))
        
        dbg("[data_corplot_data()] done!")

        res = list(Rtop=Rtop, offset=offset, klr=klr)
        return(res)
    })
    
    data_corplot.RENDER %<a-% reactive({

        dbg("[data_corplot.RENDER] reacted")
        
        res <- data_corplot_data()
        if(is.null(res)) return(NULL)
        
        par(mar=c(6,4,2,1), mgp=c(2.2,0.8,0))
        
        mar=MARGINS1
        par(mar=mar, mgp=c(1.5,0.5,0))
        
        barplot( t(res$Rtop) - res$offset, col=res$klr, border=NA, ##horiz=TRUE, 
                las=3, cex.names=0.73, ##names.arg=rep(NA,nrow(R)),
                offset = res$offset, ylab="cumulative correlation (r)"            
                ##cex.main=1.2, main="cumulative correlation\nwith other data sets"
                )
        if(!is.null(colnames(res$Rtop))) {
            legend("topright", legend=rev(colnames(res$Rtop)), fill=rev(res$klr),
                   cex=0.8, y.intersp=0.8)
        }
        dbg("[data_corplot.RENDER] done!")
        
    })
    
    ##data_corplot_module <- plotModule(
    ##id="data_corplot", ns=ns, 
    callModule(
        plotModule, "data_corplot",
        func = data_corplot.RENDER,
        func2 = data_corplot.RENDER,
        info.text = data_corplot_text,
        height = imgH, pdf.width=9, pdf.height=6,
        label="f", title="Cumulative correlation"
    )
    ##output <- attachModule(output, data_corplot_module) 
    
    ##----------------------------------------------------------------------
    ## Gene information
    ##----------------------------------------------------------------------

    data_geneInfo.RENDER  %<a-% reactive({

        require(KEGG.db)
        require(GO.db)
        require(org.Hs.eg.db)
        
        gene="A1BG-AS1"
        gene = "CD4"
        gene <- input$search_gene
        gene = toupper(sub(".*:","",gene))

        eg = "1017"
        eg = names(which(as.list(org.Hs.egSYMBOL)==gene))
        eg <- mget(gene, env=org.Hs.egSYMBOL2EG, ifnotfound=NA)[[1]]
        if(is.na(eg)) eg <- mget(gene, env=org.Hs.egALIAS2EG, ifnotfound=NA)[[1]]
        eg
        eg = eg[1]
        if(is.null(eg) || length(eg)==0) return(NULL)
        
        output = "(gene info not available)"
        if(length(eg)>0 && !is.na(eg)) {
            ##as.list(org.Hs.egSYMBOL)[[eg]]
            info <- getHSGeneInfo(eg)  ## defined in pgx-functions.R
            if(input$data_geneinfo) {
                info <- c(info, getMyGeneInfo(eg, fields="summary"))  ## defined in pgx-functions.R
            }
            output <- c()
            for(i in 1:length(info)) {
                xx <- paste(info[[i]], collapse=", ")
                output[[i]] <- paste0("<b>",names(info)[i],"</b>: ",xx)
            }
            output <- paste(output, collapse="<p>")
        }    
        ##output <- paste0("<div style='background-color: #dde6f0;'>",output,"</div>")
        ##div(HTML(output), class="gene-info-output", style="overflow: auto; height: 260px;")
        div(HTML(output), class="gene-info-output", style="overflow: auto;")
        ##div(HTML(output), class="gene-info-output")
    })

    callModule(
        plotModule, "data_geneInfo",
        plotlib = "generic",
        func = data_geneInfo.RENDER,
        func2 = data_geneInfo.RENDER,
        renderFunc = "renderUI", outputFunc = "htmlOutput",
        just.info = FALSE, no.download = TRUE,
        label="e", info.text = data_geneInfo_text,
        title = "Gene info",
        options = tagList(
            tipify( checkboxInput(ns('data_geneinfo'),'gene summary',FALSE),
                   "Provide a summary for the selected gene.", placement="top")
        ),
        height = c(260,600), width=c('auto',800)
    )

    ##----------------------------------------------------------------------
    ##                     Interface
    ##----------------------------------------------------------------------
    dataview_caption1 = "<b>Gene plots.</b> <b>(a)</b> t-SNE of samples colored by expression of selected gene. <b>(b)</b> Abundance/expression of selected gene across groups. <b>(c)</b> Top correlated genes. Level of grey corresponds to absolute expression of the gene. <b>(d)</b> Average rank of the selected gene compared to other genes. <b>(e)</b> Further information about the selected gene from public databases. <b>(f)</b> Barplot showing cumulative correlation in other datasets. <b>(g)</b> Tissue expression of selected gene."

    output$plotsUI <- renderUI({
        fillCol(
            height = fullH,
            flex = c(1,0.2,1,NA),
            fillRow( 
                flex = c(1,1,1,1), id = "data_genePlots_row1",
                height = rowH, ## width=1600, 
                ##moduleWidget(data_genePlots_tsne_module, ns=ns, height=imgH),
                ##moduleWidget(data_genePlots_barplot_module, ns=ns, height=imgH),
                ##moduleWidget(data_genePlots_correlationplot_module, ns=ns, height=imgH),
                ##moduleWidget(data_genePlots_averageRankPlot_module, ns=ns, height=imgH)
                plotWidget(ns("data_genePlots_tsne")),
                plotWidget(ns("data_genePlots_barplot")),
                plotWidget(ns("data_genePlots_correlationplot")),
                plotWidget(ns("data_genePlots_averageRankPlot"))
            ),
            br(),
            fillRow( 
                flex = c(1,1.4,1.6), id = "data_genePlots_row2",
                height = rowH, ## width=1600, 
                ## fillCol(
                ##     flex = c(NA, 0.01, 1),
                ##     data_geneInfo_buttons, br(),
                ##     htmlOutput(ns('data_geneInfo'), class="gene-info-output")
                ## ),
                ##moduleWidget(data_corplot_module, ns=ns, height=imgH),
                ##moduleWidget(data_tissueplot_module, ns=ns, height=imgH)
                plotWidget(ns("data_geneInfo")),
                plotWidget(ns("data_corplot")),
                plotWidget(ns("data_tissueplot"))
            ),
            div(HTML(dataview_caption1), class="caption")
        )
    })

    ##dragula("data_genePlots_row1")
    ##dragula("data_genePlots_row2")
    ##dragula(c("data_genePlots_row1","data_genePlots_row2"))
    
    ##----------------------------------------------------------------------
    ##                     Info messsages for Counts
    ##----------------------------------------------------------------------
    
    counts_tab_barplot_text=paste0('Barplot of the total number of counts (abundance) for each group. The samples (or cells) can be grouped/ungrouped in the ',menu_grouped, ' setting uder the main <i>Options</i>.')
    counts_tab_boxplot_text=paste0('Boxplot of the total number of counts (abundance) for each group. The samples (or cells) can be grouped/ungrouped in the ',menu_grouped, ' setting uder the main <i>Options</i>.')
    counts_tab_histplot_text=paste0('Histogram of the total number of counts (abundance) for each group. The samples (or cells) can be grouped/ungrouped in the ',menu_grouped, ' setting uder the main <i>Options</i>.')
    counts_tab_abundanceplot_text=paste0('Barplot showing the percentage of counts in terms of major gene types such as CD molecules, kinanses or RNA binding motifs for each group. The samples (or cells) can be grouped/ungrouped in the ',menu_grouped, ' setting uder the main <i>Options</i>.')
    counts_tab_average_countplot_text=paste0('Barplot showing the average count levels of major gene types such as CD molecules, kinanses or RNA binding motifs for each group. The samples (or cells) can be grouped/ungrouped in the ',menu_grouped, ' setting uder the main <i>Options</i>.')
    
    
    ##----------------------------------------------------------------------
    ##                     Count information barplot
    ##----------------------------------------------------------------------
    MARGINS2 = c(9,3.5,2,1)
    MARGINS2 = c(8,3.5,2,1)
    
    counts_tab_barplot.RENDER %<a-% reactive({
        res = getCountsTable()
        if(is.null(res)) return(NULL)
        req(input$data_groupby)
                                        #par(mar=c(20,6,10,6), mgp=c(2.2,0.8,0))
        par(mar=c(8,4,1,2), mgp=c(2.2,0.8,0))
        par(mar=MARGINS2, mgp=c(2.2,0.8,0))
        
        ylab = "counts (million)"
        if(input$data_groupby!="<ungrouped>") {
            ylab = "average group counts (million)"
        }
        subtt0=""
        if(!is.null(res$subtt)) subtt0 <- paste0("\n(",paste(res$subtt,collapse=","),")")
        ## ---- xlab ------ ###
        names.arg = names(res$total.counts)
        if( length(names.arg) > 20){ names.arg = "" }
        cex.names <- ifelse(length(names.arg)>10,0.8,0.9)
        ## ---- xlab ------ ###
        barplot(res$total.counts/1e6, las=3, border = NA,
                ##main=paste("total counts",subtt0), cex.main=1.7,
                col=rgb(0.2,0.5,0.8,0.8), #col="grey40", 
                                        #cex.names=res$cx1+0.01, cex.lab=1, ylab=ylab,
                cex.names=cex.names, cex.lab=1, ylab=ylab,
                ylim=c(0,max(res$total.counts)/1e6)*1.1,
                names.arg=names.arg)
    })

    ##counts_tab_barplot_module <- plotModule(
    ##id="counts_tab_barplot", ns=ns,
    callModule(
        plotModule, "counts_tab_barplot",
        func = counts_tab_barplot.RENDER,
        func2 = counts_tab_barplot.RENDER,        
        info.text = counts_tab_barplot_text,
        height=imgH, pdf.width=7, pdf.height=6, ## res=45,
        label="a",title='Total counts'
    )
    ## output <- attachModule(output, counts_tab_barplot_module)

    ##----------------------------------------------------------------------
    ##                     Count information boxplot
    ##----------------------------------------------------------------------

    counts_tab_boxplot.RENDER %<a-% reactive({
        res = getCountsTable()
        if(is.null(res)) return(NULL)
        ##par(mar=c(3,3,3,3), mgp=c(2.4,0.7,0), oma=c(1,1,1,1)*0.2 )   
        par(mar=c(8,4,1,2), mgp=c(2.2,0.8,0))
        par(mar=MARGINS2, mgp=c(2.2,0.8,0))
        
        ## ---- xlab ------ ###
        xaxt="l"
        names.arg = colnames(res$log2counts[res$jj,])
        if( length(names.arg) > 20){ names.arg = rep("",length(names.arg)); xaxt="n"}
        cex.names <- ifelse(length(names.arg)>10,0.8,0.9)
        boxplot(res$log2counts[res$jj,], col=rgb(0.2,0.5,0.8,0.4), #col=rgb(0.2,0.5,0.8,0.3), #col="grey70", 
                                        #main="counts distribution", cex.main=1.6,
                ## cex.names=res$cx1+0.1,
                names = names.arg, cex.axis=cex.names,#border=rgb(0.2,0.5,0.8,0.8),
                border = 	rgb(0.824,0.824,0.824,0.9),xaxt=xaxt, 
                las=3, cex.lab=1, ylab="counts (log2)", outline=FALSE, varwidth = FALSE)
    })

    ##counts_tab_boxplot_module <- plotModule(
    ##    id="counts_tab_boxplot", ns=ns,
    callModule(
        plotModule, "counts_tab_boxplot",
        func = counts_tab_boxplot.RENDER,
        func2 = counts_tab_barplot.RENDER,
        info.text = counts_tab_boxplot_text,
        height=imgH, pdf.width=7, pdf.height=6, ## res=50,
        label="b", title='Counts distribution'
    )
    ##output <- attachModule(output, counts_tab_boxplot_module)

    ##----------------------------------------------------------------------
    ##                     Count information histogram
    ##----------------------------------------------------------------------

    counts_tab_histplot.RENDER %<a-% reactive({
        res = getCountsTable()
        if(is.null(res)) return(NULL)
                                        #par(mfrow=c(2,3), mar=c(9,4,3,1.5), mgp=c(2.4,0.7,0), oma=c(1,1,1,1)*0.2 )  
        par(mar=c(8,4,1,2), mgp=c(2.2,0.8,0))
        par(mar=MARGINS2, mgp=c(2.2,0.8,0))
        
        n=1000
        gx.hist <- function(gx, n=1000, main="",ylim=NULL) {
            jj <- 1:nrow(gx)
            if(length(jj)>n) jj <- sample(jj,n,replace=TRUE)
            h0 <- hist(as.vector(c(gx[jj],min(gx),max(gx))),
                       breaks=120, main=main, border=FALSE,
                                        #col=rgb(0.2,0.5,0.8,0.5), 
                       col="grey", 
                       freq=FALSE, ## ylim=ylim,
                       xlim=c(min(gx),max(gx)),
                       xlab="expression (log2)", 
                       ##cex.names=res$cx1+0.1,
                       cex.lab=1)
            i = 1
            for(i in 1:ncol(gx)) {
                h1 <- hist(gx[jj,i], breaks=h0$breaks,plot=FALSE)
                ##lines( h0$mids, h1$density, col=i+1 )
                lines( h0$mids, h1$density, col="black", lwd=0.5 )
            }
        }
        gx.hist(gx=res$log2counts[res$jj,], n=2000) #, main="histogram")
    })

    ##counts_tab_histplot_module <- plotModule(
    ##id="counts_tab_histplot", ns=ns,
    callModule(
        plotModule, "counts_tab_histplot",
        func = counts_tab_histplot.RENDER,
        func2 = counts_tab_histplot.RENDER,
        info.text = counts_tab_histplot_text,
        height=imgH, pdf.width=7, pdf.height=6, ## res=50,
        label="c", title='Counts histogram'
    )
    ##output <- attachModule(output, counts_tab_histplot_module)

    ##----------------------------------------------------------------------
    ##  Count information abundance of major gene types
    ##----------------------------------------------------------------------

    counts_tab_abundanceplot.RENDER %<a-% reactive({
        res = getCountsTable()
        if(is.null(res)) return(NULL)
                                        #par(mar=c(6,4,0,4), mgp=c(2.2,0.8,0))
        par(mar=c(6,4,0,4))
        par(mar=c(7,4,0,2))
        par(mar=MARGINS2, mgp=c(2.2,0.8,0))
        
        klr <- colorRampPalette(c(rgb(0.2,0.5,0.8,0.8), rgb(0.2,0.5,0.8,0.1)), alpha = TRUE)(nrow(res$prop.counts))
        
        ## ----------- gene type counts
        ymax = max(colSums(res$prop.counts, na.rm=TRUE))
        ## ---- xlab ------ ###
        names.arg = colnames(res$prop.counts)
        if( length(names.arg) > 20){ names.arg = rep("",length(names.arg)) }
        cex.names <- ifelse(length(names.arg)>10,0.8,0.9)
        ## ---- xlab ------ ###
        barplot(res$prop.counts, las=3, #main="abundance of major gene types", cex.main=1.6,
                                        #cex.names=res$cx1+0.04, 
                cex.lab=1.0, border = NA,
                ylim=c(0,ymax)*1.6, ylab="abundance (%)",
                names.arg=names.arg, cex.names=cex.names,
                col = klr)
                                        # legend("topleft", legend=rev(rownames(res$prop.counts)),
                                        #        fill=rev(grey.colors(nrow(res$prop.counts))),
                                        #        cex=1, y.intersp=0.4, bty="n")
        leg <- legend("topleft", legend=rev(rownames(res$prop.counts)),
                                        #fill=rev(grey.colors(nrow(res$prop.counts))),
                      fill=rev(klr),cex=1, y.intersp=0.75, bty="n", plot = FALSE)
        leftx <- leg$rect$left*0.9
        rightx <- leg$rect$right*0.9
        topy <- leg$rect$top
        bottomy <- leg$rect$bottom
        legend(x = c(leftx, rightx), y = c(topy, bottomy), 
               legend=rev(rownames(res$prop.counts)),
                                        #fill=rev(grey.colors(nrow(res$prop.counts))),
               fill=rev(klr), bty="n", cex=0.9, y.intersp=0.75)

    })

    ##counts_tab_abundanceplot_module <- plotModule(
    ##id="counts_tab_abundanceplot", ns=ns,
    callModule(
        plotModule, "counts_tab_abundanceplot",        
        func = counts_tab_abundanceplot.RENDER,
        func2 = counts_tab_abundanceplot.RENDER,
        info.text = counts_tab_abundanceplot_text,
        height=imgH, pdf.width=10, pdf.height=6, ## res=50,
        label="d",title='Abundance of major gene types'
    )
    ##output <- attachModule(output, counts_tab_abundanceplot_module) 

    ##----------------------------------------------------------------------
    ## Count information average count by gene type
    ##----------------------------------------------------------------------
    counts_tab_average_countplot.RENDER %<a-% reactive({
        res = getCountsTable()
        if(is.null(res)) return(NULL)
                                        #par(mar=c(6,4,0,4), mgp=c(2.2,0.8,0))
        par(mar=c(6,4,0,4))
        par(mar=c(7,4,0,2))
        par(mar=MARGINS2, mgp=c(2.2,0.8,0))
        
        klr <- colorRampPalette(c(rgb(0.2,0.5,0.8,0.8), rgb(0.2,0.5,0.8,0.1)), alpha = TRUE)(nrow(res$avg.counts))
                                        # klr <- grey.colors(nrow(res$avg.counts))
        
        ## ---- xlab ------ ###
        names.arg = colnames(res$avg.counts)
        if( length(names.arg) > 20){ names.arg = rep("",length(names.arg)) }
        cex.names <- ifelse(length(names.arg)>10,0.8,0.9)
        ## ---- xlab ------ ###
        ymax = max(colSums(res$avg.counts, na.rm=TRUE))
        barplot(res$avg.counts, las=3, #main="average counts by gene type", cex.main=1.6,
                                        #cex.names=res$cx1+0.04, 
                border=NA, cex.lab=1.0, 
                names.arg=names.arg, cex.names=cex.names,
                ylim=c(0,ymax)*1.6, ylab="average count", col=klr)
        leg <- legend("topleft", legend=rev(rownames(res$avg.counts)),
                                        #fill=rev(grey.colors(nrow(res$avg.counts))),
                      fill=rev(klr), cex=1, y.intersp=0.75, bty="n", plot = FALSE)
        leftx <- leg$rect$left*0.9
        rightx <- leg$rect$right*0.9
        topy <- leg$rect$top
        bottomy <- leg$rect$bottom
        legend(x = c(leftx, rightx), y = c(topy, bottomy),
               legend=rev(rownames(res$avg.counts)),
                                        #fill=rev(grey.colors(nrow(res$avg.counts))),
               fill=rev(klr), cex=0.9, y.intersp=0.75, bty="n")
    })

    ##counts_tab_average_countplot_module <- plotModule(
    ##id="counts_tab_average_countplot", ns=ns,
    callModule(
        plotModule, "counts_tab_average_countplot",
        func = counts_tab_average_countplot.RENDER,
        func2 = counts_tab_average_countplot.RENDER,
        info.text = counts_tab_average_countplot_text,
        height=imgH, pdf.width=10, pdf.height=6, ##res=50,
        label="e",title='Average count by gene type'
    )
    ##output <- attachModule(output, counts_tab_average_countplot_module) 

    getCountsTable <- reactive({
        ngs = inputData()
        req(ngs)
        
        validate(need("counts" %in% names(ngs), "no 'counts' in object."))    
        subtt=NULL

        samples = colnames(ngs$X)
        samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        nsamples = length(samples)
        if("counts" %in% names(ngs)) {
            counts = ngs$counts[,samples,drop=FALSE]
        } else {
            cat("WARNING:: no counts table. estimating from X\n")
            counts = pmax(2**ngs$X - 1,0)
            k = grep("lib.size",colnames(ngs$samples))[1]
            if(length(k)>0) {
                libsize = ngs$samples[colnames(counts),k]
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
        gr = ngs$Y[samples,grpvar]
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
        gg = ngs$genes[rownames(counts),]$gene_name
        tt = ngs$genes[rownames(counts),]$gene_title
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

        res = list(total.counts=total.counts,subtt=subtt, cx1=cx1,
                   log2counts=log2counts, jj=jj, prop.counts=prop.counts,
                   avg.counts=avg.counts)

    })

    dataview_counts_caption = "<b>Counts distribution</b>. Plots associated with the counts, abundance or expression levels across the samples/groups.  <b>(a)</b> Total counts per sample or average per group.  <b>(b)</b> Distribution of total counts per sample/group. The center horizontal bar correspond to the median.  <b>(c)</b> Histograms of total counts distribution per sample/group. <b>(d)</b> Abundance of major gene types per sample/group. <b>(e)</b> Average count by gene type per sample/group."

    output$countsUI <- renderUI({
        fillCol(
            flex = c(1,1,NA),
            height = fullH,
            fillRow(
                flex = c(1,1,1), id = "counts_tab_row1", height=rowH,
                ##moduleWidget(counts_tab_barplot_module, ns=ns, height=imgH),
                ##moduleWidget(counts_tab_boxplot_module, ns=ns, height=imgH),
                ##moduleWidget(counts_tab_histplot_module, ns=ns, height=imgH)
                plotWidget(ns("counts_tab_barplot")),
                plotWidget(ns("counts_tab_boxplot")),
                plotWidget(ns("counts_tab_histplot"))
            ),
            fillRow(
                flex = c(1,1), id = "counts_tab_row2", height=rowH,
                ##moduleWidget(counts_tab_abundanceplot_module, ns=ns, height=imgH),
                ##moduleWidget(counts_tab_average_countplot_module, ns=ns, height=imgH)
                plotWidget(ns("counts_tab_abundanceplot")),
                plotWidget(ns("counts_tab_average_countplot"))
            ),
            div(HTML(dataview_counts_caption), class="caption")
        )
    })
    ##dragula(c("counts_tab_row1","counts_tab_row2"))
    
    ##================================================================================
    ##======================  Raw counts/abundance table =============================
    ##================================================================================

    dropdown_search_gene='<code>Search gene</code>'
    menu_grouped='<code>grouped</code>'
    menu_options='<code>Options</code>'
    
    data_rawdataTable_text = paste0('Under the <strong>gene table </strong>, the average expression values of genes across the groups can be read. The samples (or cells) can be ungrouped by unclicking the ',menu_grouped, ' in the main <i>Options</i> to see the exact expression values per sample (or cell).', 'The genes in the table are ordered by the correlation (<b>rho</b> column) with respect to the gene selected by users from the ',dropdown_search_gene, ' setting. <b>SD</b> column reports the standard deviation of expression across samples (or cells).')
    
    data_rawdataTable.RENDER <- reactive({
        ## get current view of raw_counts
        ngs = inputData()
        req(ngs)
        req(input$data_groupby)

        dbg("[data_rawdataTable.RENDER] reacted")
        
        pp <- rownames(ngs$X)
        if(input$data_type=="counts") {
            x <- ngs$counts[pp,]
        } else if(input$data_type=="CPM") {
            ##x <- 2**ngs$X
            x <- edgeR::cpm(ngs$counts[pp,])
        } else {
            ## log2CPM
            x <- ngs$X[pp,]
        }
        x0=x
        
        ##------------------ select samples
        x <- x[,]
        samples <- colnames(ngs$X)
        samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        samples <- intersect(colnames(x),samples)
        x <- x[,samples,drop=FALSE]
        
        gene = "CD3E"
        gene = "CCR6"
        gene = input$search_gene
        if(is.null(gene) | gene=="" | is.na(gene)) return(NULL)
        ##xgene = sub(".*:","",rownames(x))

        ## Quickly (?) calculated correlation to selected gene
        rho = NULL
        if(1) {
            k=1
            ##k <- which(xgene==gene)
            ##logx <- log2(1 + edgeR::cpm(ngs$counts))
            ##logx <- logx[rownames(x),]
            logx <- ngs$X[rownames(x),]
            xgenes <- ngs$genes[rownames(x),"gene_name"]
            k <- which(xgenes==gene)
            rho = cor( t(logx[,samples]), logx[k,samples], use="pairwise")[,1]
            rho = round(rho[rownames(x)], digits=3)
            sdx = round(apply(logx[,samples],1,sd),digits=3)
        }
        
        ##if(input$data_sampling=="grouped") {
        ##do.grouped <- input$data_grouped
        do.grouped <- (input$data_groupby != "<ungrouped>")
        if(length(samples)>500) do.grouped <- TRUE
        if(do.grouped) {
            grpvar <- input$data_groupby
            group = ngs$Y[colnames(x),grpvar]
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
        ##rownames(x) = sub(".*:","",rownames(x))
        xgenes <- ngs$genes[rownames(x),"gene_name"]
        gene.title <- GENE.TITLE[toupper(xgenes)]
        gene.title <- substring(gene.title,1,50)
        x = data.frame( gene=xgenes, title=gene.title, rho=rho, SD=sdx,
                       as.matrix(x), check.names=FALSE)
        if(!is.null(rho)) {
            x = x[order(-abs(x$rho)),,drop=FALSE]
        } else {
            x = x[order(-x$SD),,drop=FALSE]
        }

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
        
        dbg("[data_rawdataTable.RENDER] rendering N=",nrow(x),"rows...")
        
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
                DT::formatStyle(colnames(x),
                                background = styleColorBar(c(0,x99), 'lightblue'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')
    })

    data_rawdataTable_caption = "<b>Gene table.</b> The table shows the gene expression values per sample, or average expression values across the groups. The column 'rho' reports the correlation with the gene selected in 'Search gene' in the left side bar."

    ## data_rawdataTable_module <- tableModule(
    ##     id = "data_rawdataTable", ns=ns,
    ##     func = data_rawdataTable.RENDER,
    ##     title = "Gene expression table",
    ##     info.text = data_rawdataTable_text,
    ##     caption = data_rawdataTable_caption
    ## )
    ## output <- attachModule(output, data_rawdataTable_module) 
    data_rawdataTable <- callModule(
        tableModule, "data_rawdataTable",
        func = data_rawdataTable.RENDER,
        title = "Gene expression table",
        info.text = data_rawdataTable_text,
        caption = data_rawdataTable_caption
    )

    output$genetableUI <- renderUI({
        fillCol(
            flex = c(1),
            height = fullH,
            ##moduleWidget(data_rawdataTable_module, outputFunc="dataTableOutput",ns=ns)
            tableWidget(ns("data_rawdataTable"))
        )
    })
    
    ##================================================================================
    ##================================= Samples ======================================
    ##================================================================================

    data_phenoClustering.RENDER <- reactive({
        ngs = inputData()
        req(ngs)
        dbg("[data_phenoClustering.RENDER] reacted")

        annot <- ngs$samples
        samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        annot <- annot[samples,,drop=FALSE]
        annot.ht <- ifelse( ncol(annot) > 10, 3, 5)

        do.clust <- input$data_phenoclustsamples
        plt <- pgx.plotPhenotypeMatrix0(
            annot, annot.ht=annot.ht, cluster.samples=do.clust)
        ## plt <- plt %>% config(displayModeBar = FALSE)
        plt
    })

    data_phenoClustering_opts <- tagList(
        tipify( checkboxInput(ns('data_phenoclustsamples'),'cluster samples',TRUE),
               "Cluster samples.", placement="top")        
    )
        
    data_phenoClustering_caption = "<b>Phenotype clustering.</b> Clustered heatmap of sample information (phenotype data)."
    callModule(
        plotModule,
        id = "data_phenoClustering", label="a",
        func = data_phenoClustering.RENDER,
        func2 = data_phenoClustering.RENDER,        
        ## plotlib = "iheatmapr", 
        title = "Phenotype clustering",
        ##info.text = "Sample information table with information about phenotype of samples.",
        options = data_phenoClustering_opts,
        height = c(320,600), width = c('auto',1200), res=c(68,75),
        pdf.width=10, pdf.height=6 
    )
    ##output <- attachModule(output, data_sampleTable_module) 
    
    data_sampleTable.RENDER <- reactive({
        ## get current view of raw_counts
        ngs = inputData()
        req(ngs)

        dbg("[data_sampleTable.RENDER] reacted")
        
        ##if(is.null(input$data_samplefilter)) return(NULL)    
        dt <- NULL
        samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        dt <- ngs$samples[samples,,drop=FALSE]    
        DT::datatable( dt,
                      class = 'compact cell-border stripe hover',
                      rownames = TRUE,
                      extensions = c('Buttons','Scroller'),
                      selection = list(mode='single', target='row', selected=1),
                      options=list(
                          dom = 'lfrtip', 
                          ##pageLength = 60, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scroller=TRUE, scrollX = TRUE,
                          scrollY = 0.35*tabH,
                          deferRender=TRUE
                      )) %>%
            DT::formatStyle(0, target='row', fontSize='12px', lineHeight='70%') 
        
    })

    data_sampleTable_caption="<b>Sample information.</b> This table contains available phenotype information about the samples. Phenotype variables starting with a 'dot' (e.g. '.cell cycle' and '.gender' ) have been estimated from the data."
    ##data_sampleTable_module <- tableModule(
    ##id="data_sampleTable", ns=ns,
    data_sampleTable <- callModule(
        tableModule, "data_sampleTable", label="b",
        func = data_sampleTable.RENDER,
        title = "Sample information",
        info.text = "Sample information table with information about phenotype of samples.",
        ##options = data_sampleTable_opts,
        ## caption = data_sampleTable_caption
    )
    ##output <- attachModule(output, data_sampleTable_module) 

    sampletableUI_caption <- paste(
        ## "<b>Phenotype clustering and sample information table.</b>",
        "<b>(a)</b>",data_phenoClustering_caption,
        "<b>(b)</b>",data_sampleTable_caption
    )
    output$sampletableUI <- renderUI({
        fillCol(
            flex = c(1,1,NA),
            height = fullH,
            ##moduleWidget(data_sampleTable_module, outputFunc="dataTableOutput", ns=ns)
            div( plotWidget(ns("data_phenoClustering")), style="overflow: auto;"),
            tableWidget(ns("data_sampleTable")),
            div(HTML(sampletableUI_caption), class="caption")
        )
    })
        
    ##================================================================================
    ##================================= CONTRASTS ====================================
    ##================================================================================

    data_contrastTable.RENDER <- reactive({
        ## get current view of raw_counts
        ngs = inputData()
        req(ngs)

        dbg("[data_contrastTable.RENDER] reacted")
        
        ##if(is.null(input$data_samplefilter)) return(NULL)    
        dt <- NULL
        samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        names(ngs$model.parameters)
        if(input$data_ctbygroup=="group") {
            ct <- ngs$model.parameters$contr.matrix
            kk <- which(rownames(ct) %in% ngs$samples[samples,"group"])
            dt <- ct[kk,,drop=FALSE]
        } else {
            dt <- ngs$model.parameters$exp.matrix[samples,,drop=FALSE]
        }
        dt <- sign(dt)
        colnames(dt) <- sub("[_. ]vs[_. ]","\nvs ",colnames(dt))
        dt[dt==0] <- NA

        ## looks better this way
        dt1 <- dt
        if(ncol(dt)<8) {
            dt1 <- cbind(dt,NA,NA,NA,NA,NA,NA,NA,NA)[,1:8]
        }
        
        DT::datatable( dt1,
                      class = 'compact cell-border stripe hover',
                      rownames = TRUE,
                      extensions = c('Buttons','Scroller'),
                      selection = list(mode='single', target='row', selected=1),
                      options=list(
                          dom = 'lfrtip', 
                          ##pageLength = 60, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scroller=TRUE, scrollX = TRUE, scrollY = tabH,
                          deferRender=TRUE
                      )) %>%
            DT::formatStyle(0, target='row', fontSize='12px', lineHeight='70%') %>%
                DT::formatStyle(colnames(dt),
                                ##background = styleColorBar(c(0,3), 'lightblue'),
                                background = color_from_middle(c(-1,1), 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')
        
        
    })

    data_contrastTable_info = "Table summarizing the contrasts of all comparisons. Here, you can check which samples belong to which groups for the different comparisons. Non-zero entries '+1' and '-1' correspond to the group of interest and control group, respectively. Zero or empty entries denote samples not use for that comparison."
    
    data_contrastTable_caption = "<b>Contrast table</b> summarizing the contrasts of all comparisons. Non-zero entries '+1' and '-1' correspond to the group of interest and control group, respectively. Zero or empty entries denote samples not use for that comparison."

    data_contrastTable_opts = tagList(
        tipify( radioButtons(ns('data_ctbygroup'), "Show by:", choices=c("group","sample")),
               "Show contrasts by group or by samples.",
               placement="right", options = list(container = "body"))
    )
    
    ##data_contrastTable_module <- tableModule(
    ##id = "data_contrastTable", ns=ns,
    data_contrastTable <- callModule(
        tableModule, "data_contrastTable",
        func = data_contrastTable.RENDER,
        options = data_contrastTable_opts,
        title = "Contrast table",
        info.text = data_contrastTable_info,
        caption = data_contrastTable_caption
    )
    ##output <- attachModule(output, data_contrastTable_module) 

    output$contrasttableUI <- renderUI({
        fillCol(
            flex = c(1), height = fullH,
            ##moduleWidget(data_contrastTable_module, outputFunc="dataTableOutput", ns=ns)
            tableWidget(ns("data_contrastTable"))
        )
    })
    
    ##================================================================================
    ## Resource info (dev)
    ##================================================================================
    
    datatable_timings.RENDER <- reactive({
        ngs <- inputData()
        req(ngs)

        dbg("[datatable_timings.RENDER] reacted")
        
        ##if(is.null(ngs$timings)) return(NULL)
        D <- data.frame()
        if(!is.null(ngs$timings)) {
            D <- round(ngs$timings[,1:3],digits=3)
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
    ##datatable_timings_module <- tableModule(
    ##id="datatable_timings", ns=ns,
    datatable_timings <- callModule(
        tableModule, "datatable_timings",
        func=datatable_timings.RENDER,
        info.text = datatable_timings_text,
        options = NULL, title='Timings'
    )
    ## output <- attachModule(output, datatable_timings_module)
    
    datatable_objectdims.RENDER <- reactive({
        ngs <- inputData()
        req(ngs)
        
        dims1 <- lapply( ngs, dim)
        lens <- sapply( ngs, length)
        dims2 <- t(sapply( ngs[which(!sapply(dims1,is.null)) ], dim))
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

    ##datatable_objectdims_module <- tableModule(
    ##    id="datatable_objectdims", ns=ns,
    datatable_objectdims <- callModule(
        tableModule, "datatable_objectdims",        
        func = datatable_objectdims.RENDER,
        info.text = datatable_objectdims_text,
        options = NULL, title='Object dimensions'
    )
    ##output <- attachModule(output, datatable_objectdims_module)    
    
    datatable_objectsize.RENDER <- reactive({
        ngs <- inputData()
        req(ngs)
        objsize <- sapply(ngs,object.size)
        objsize <- round( objsize/1e6, digits=2)
        D = data.frame( object=names(ngs), "size.Mb"=objsize, check.names=FALSE)
        DT::datatable( D, rownames=FALSE,
                      options = list(dom='t', pageLength = 50), 
                      class = 'compact cell-border stripe hover') %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')
    })
    
    datatable_objectsize_text = "This table provides information about  about the memory sizes of objects"

    ##datatable_objectsize_module <- tableModule(
    ##    id="datatable_objectsize", ns=ns,
    datatable_objectsize <- callModule(
        tableModule, "datatable_objectsize",            
        func = datatable_objectsize.RENDER,
        options = NULL, title='Object sizes',
        info.text = datatable_objectsize_text
        ## caption = datatable_objectsize_caption
    )
    ## output <- attachModule(output, datatable_objectsize_module)

    resourceinfo_caption="<b>Resource information.</b> Details about the execution times of the methods, dimensions and memory sizes of objects."
    
    output$resourceinfoUI <- renderUI({    
        fillCol(
            flex = c(1,NA),
            height = fullH,
            fillRow(
                flex = c(5,1, 2,1, 1.5, 2), ## width = 600,
                ##moduleWidget(datatable_timings_module, outputFunc="dataTableOutput", ns=ns),
                tableWidget(ns("datatable_timings")),                
                br(),
                ##moduleWidget(datatable_objectdims_module, outputFunc="dataTableOutput", ns=ns),
                tableWidget(ns("datatable_objectdims")),                
                br(),
                ##moduleWidget(datatable_objectsize_module, outputFunc="dataTableOutput", ns=ns),
                tableWidget(ns("datatable_objectsize")),                
                br()
            ),
            div(HTML(resourceinfo_caption),class="caption")
        )
    })


}
