SignatureInputs <- function(id) {
    ns <- NS(id)  ## namespace
    tagList(
        uiOutput(ns("description")),
        uiOutput(ns("inputsUI"))
    )
}

SignatureUI <- function(id) {
    ns <- NS(id)  ## namespace
    fillRow(
        flex = c(1.5,0.05,1),
        height = 750,
        tabsetPanel(
            tabPanel("Enrichment",uiOutput(ns("sig_enplots_UI"))),
            tabPanel("Overlap/similarity",uiOutput(ns("sig_overlapAnalysis_UI"))),
            tabPanel("Markers",uiOutput(ns("sig_markers_UI")))
        ),
        br(),
        tabsetPanel(
            tabPanel("Enrichment table",uiOutput(ns("sig_enrichmentTables_UI")))
        )
    )
}

SignatureModule <- function(input, output, session, inputData)
{
    ns <- session$ns ## NAMESPACE

    description = "<b>Signature Analysis.</b> Users can test their gene signature by
calculating an enrichment score. Upload your own gene list, or select
a contrast which then takes the top differentially expressed genes as
signature."
    output$description <- renderUI(HTML(description))


sig_infotext =
    "In the <strong>Signature Analysis module</strong>, users can test their gene signature by calculating an enrichment score. They can use a sample list provided on the platform or upload their own gene list. Instead of a short list, a profile can also be selected, which is a complete gene list resulted from one of the contrasts in the analysis.

<br><br>After uploading a gene list, the <strong>Markers</strong> section produces a t-SNE plot of samples for each gene, where the samples are colored with respect to the upregulation (in red) or downregulation (in blue) of that particular gene.

<br><br>The <strong>Enrichment tab</strong> performs the enrichment analysis of the gene list against all contrasts by running the GSEA algorithm and plots enrichment outputs. The enrichment statistics can be found in the corresponding table

<br><br>Under the <strong>Overlap/similarity tab</strong>, users can find the similarity of their gene list with all the gene sets and pathways in the platform, including statistics such as the total number of genes in the gene set (K), the number of intersecting genes between the list and the gene set (k), the overlapping ratio of k/K, as well as the p and q values by the Fisher’s test for the overlap test.

<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=7' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>
"

    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================

    IMMCHECK.GENES = "ADORA2A ARHGEF5 BTLA CD160 CD244 CD27 CD274 CD276 CD47 CD80 CEACAM1 CTLA4 GEM HAVCR2 ICOS IDO1 LAG3 PDCD1 TNFSF4 VISTA VTCN1 TIGIT PVR CD28 CD40 CD40LG ICOSLG TNFRSF9 TNFSF9 CD70 TNFRSF4 TNFRSF18 TNFSF18 SIRPA LGALS9 ARG1 CD86 IDO2 PDCD1LG2 KIR2DL3"
    APOPTOSIS.GENES = "BAD CRADD AGT FAS BCL2 PPIF S100A9 S100A8 BBC3 BCL2L11 FADD CTSH MLLT11 TRAF7 BCL2L1 HTRA2 BNIP3 BAK1 PMAIP1 LGALS9 BID"
    CELL.CYCLE.GENES = "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP
5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA"
    style0 = "font-size: 0.9em; color: #24A; background-color: #dde6f0; border-style: none; padding:0"
    
    output$inputsUI <- renderUI({
        ui <- tagList(
            tipify( actionLink(ns("sig_info"), "Info", icon = icon("info-circle")),
                   "Show more information about this module"),
            hr(), br(), 
            conditionalPanel(
                ##condition="input.main_usermode=='PRO'",
                condition="output.main_usermode=='PRO' || output.main_usermode=='DEV'",
                tagList(
                    tipify(selectInput(ns("sig_type"),label="Signature type:",choices=c("<custom>","contrast","hallmark","KEGG")),
                           "Specify the type of signature of an interest. Users can choose between custom signature, a contrast profile, or some predefined gene sets including Hallmark and KEGG pathways.", placement="top", options = list(container = "body")),
                    conditionalPanel(
                        "input.sig_type != '<custom>'",
                        tipify(selectizeInput(ns("sig_feature"),"Signature:", choices="<custom>",
                                              selected="<custom>"),
                               "Select a specific signature group.", placement="top",
                               options = list(container = "body"))
                    ))
            ),
            br(),br(),
            tipify(textAreaInput(ns("sig_genelistUP"), "Genes:", value = IMMCHECK.GENES,
                                 rows=10, placeholder="Paste your gene list"),
                   "Paste a list of signature genes.", placement="top", options = list(container = "body")),
            ## textAreaInput("sig_genelistDN", "Signature (down):", rows=6, placeholder="Paste your gene list")
            tipify(actionButton(ns("sig_example2"),"[apoptosis] ", style=style0),
                   "Use the list of genes involved in apoptosis as a signature."),
            tipify(actionButton(ns("sig_example1"),"[immune_chkpt] ", style=style0),
                   "Use the list of genes involved in immune checkpoint as a signature."),
            tipify(actionButton(ns("sig_example3"),"[cell_cycle] ", style=style0),
                   "Use the list of genes involved in cell cycle as a signature."),
            br(),br(),
            conditionalPanel(
                ##condition="input.main_usermode=='PRO'",
                condition="output.main_usermode=='PRO' || output.main_usermode=='DEV'",
                tagList(
                    tipify( selectInput(ns('cmp_querydataset'),"Query dataset:", choices=NULL, multiple=FALSE),
                           "The query dataset to which the enrichment test should be applied. Enrichment of the selected signature will be calculated to all available contrast profiless in this query dataset.", 
                           placement="top", options = list(container = "body"))
                )
            )
        )

        if(DEV.VERSION) {
            uix <- tagList(
                br(),br(),hr(),h6("Developer options:"),
                radioButtons('sig_ssstats','ss-stats:',c("rho","gsva","grp.gsva","rho+gsva","rho+grp.gsva"),
                             inline=TRUE)
            )
            ui <- c(ui, uix)
        }
        ui
    })

    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    observeEvent( input$sig_info, {
        showModal(modalDialog(
            title = HTML("<strong>Signature Analysis Module</strong>"),
            HTML(sig_infotext),
            easyClose = TRUE, size="l"))
    })

    ##------------------------ observe/reactive function  -----------------------------
    require(visNetwork)
    require(igraph)
    require(parallel)
    require(fgsea)
    require(GSVA)

    observeEvent(input$sig_example1, { 
        updateTextAreaInput(session,"sig_genelistUP", value=IMMCHECK.GENES)
    })
    observeEvent(input$sig_example2, { 
        updateTextAreaInput(session,"sig_genelistUP", value=APOPTOSIS.GENES)
    })
    observeEvent(input$sig_example3, { 
        updateTextAreaInput(session,"sig_genelistUP", value=CELL.CYCLE.GENES)
    })

    observe({
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)
        sig_type="contrast"
        sig_type <- input$sig_type
        if(is.null(sig_type)) sig_type <- "<custom>"

        if(sig_type=="contrast") {
            contr <- sort(names(ngs$gx.meta$meta))
            updateSelectizeInput(session, "sig_feature", choices=contr, selected=contr[1])
        } else if(sig_type=="hallmark") {
            ## collection
            gsets <- sort(grep("HALLMARK",names(GSETS),value=TRUE))
            updateSelectizeInput(session, "sig_feature", choices=gsets, selected=gsets[1])
        } else if(sig_type=="KEGG") {
            ## collection
            gsets <- sort(grep("KEGG",names(GSETS),value=TRUE))
            updateSelectizeInput(session, "sig_feature", choices=gsets, selected=gsets[1])
        } else if(sig_type=="geneset") {
            ## all genesets... this is a bit too much for selectInput (DO NOT USE!!)
            gsets <- sort(names(GSETS))
            updateSelectizeInput(session, "sig_feature", choices=gsets, selected=gsets[1])
        } else {
            ## custom
            updateSelectizeInput(session, "sig_feature", choices="<custom>", selected="<custom>")
        }
    })


    ##================================================================================
    ##======================= REACTIVE FUNCTIONS =====================================
    ##================================================================================

    input_sig_genelistUP <- reactive({
        gg <- input$sig_genelistUP
        if(is.null(gg)) return(NULL)
        gg <- strsplit(as.character(gg), split="[, \n\t]")[[1]]
        if(length(gg)==1 && gg[1]!="") gg <- c(gg,gg)  ## hack to allow single gene....
        return(gg)
    }) %>% debounce(1000)

    
    getCurrentMarkers <- reactive({
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)

        sig_type="<custom>"
        sig_type="contrast"
        sig_type <- input$sig_type
        ##if(is.null(sig_type)) return(NULL)
        ##if(is.null(input$sig_contrast)) return(NULL)
        ##if(is.null(input$sig_feature)) return(NULL)
        req(input$sig_type, input$sig_feature)
        
        dbg("<signature:getCurrentMarkers> called\n")
        
        level = "gene"
        features = toupper(ngs$genes$gene_name)
        xfeatures = toupper(ngs$genes[rownames(ngs$X),"gene_name"])
        gset <- NULL
        if(input$sig_feature=="<custom>") {
            gset <- input_sig_genelistUP()
            if(is.null(gset) || length(gset)==0 || gset[1]=="") return(NULL)
            ##gset <- toupper(gset)        
            if(length(gset)==1) {
                gene <- sub("^[@#]","",gset[1])
                if(grepl("^@",gset[1]) && gene %in% xfeatures) {
                    ## most correlated with this genes
                    jj <- match(gene, xfeatures)  ## single gene
                    rho <- cor(t(ngs$X), ngs$X[jj,])[,1]
                    gset <- head(names(sort(abs(rho),decreasing=TRUE)),36)  ## how many?
                } else {
                    ## grep-like match
                    rx <- toupper(gset[1])
                    rx <- grep(rx, xfeatures, value=TRUE, ignore.case=TRUE)
                    gset <- rownames(ngs$X)[which(xfeatures %in% rx)]  ## all probes matching gene
                }
            }
        } else if(sig_type=="contrast" &&
                  input$sig_feature %in% names(ngs$gx.meta$meta) ) {
            contr=1
            contr <- input$sig_feature
            fx <- ngs$gx.meta$meta[[contr]]$meta.fx
            probes <- rownames(ngs$gx.meta$meta[[contr]])
            genes <- toupper(ngs$genes[probes,"gene_name"])
            top.genes <- genes[order(-fx)]
            top.genes <- intersect( top.genes, rownames(PROFILES$FC))
            top.genes <- head(top.genes,100)
            top.genes0 <- paste(top.genes,collapse=" ")
            updateTextAreaInput(session,"sig_genelistUP", value=top.genes0)
            gset <- top.genes
        } else if(input$sig_feature %in% names(GSETS)) {
            gset <- toupper(GSETS[[input$sig_feature]])
            gset <- intersect(gset, rownames(PROFILES$FC))
            gset0 <- paste(gset, collapse=" ")
            updateTextAreaInput(session,"sig_genelistUP", value=gset0)
        } else {
            return(NULL)
        }
        
        return(gset)
    })
    
    getSingleSampleEnrichment <- reactive({
        ##
        ## Calls calcSingleSampleValues() and calculates single-sample
        ## enrichment values for complete data matrix and reduced data by
        ## group (for currentmarkers)
        ##
        ##
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)
        
        dbg("getSingleSampleEnrichment:: called\n")

        ## select samples
        X = ngs$X
        sel = colnames(X)
        X <- X[,sel]
        
        ## get the signature
        gset <- strsplit(IMMCHECK.GENES,split=" ")[[1]]
        gset <- getCurrentMarkers()
        if(is.null(gset)) return(NULL)

        ##y = 1*(rownames(X) %in% gset)
        ##rownames(X)=toupper(rownames(X)); gset=toupper(gset)
        xgene <- ngs$genes[rownames(X),"gene_name"]
        y = 1*(toupper(xgene) %in% toupper(gset))
        names(y) <- rownames(X)
        table(y)
        
        dbg("getSingleSampleEnrichment:: 1")
        
        ## expression by group
        grp = ngs$samples[colnames(X),"group"]
        groups = unique(grp)
        gX <- sapply( groups, function(g) rowMeans(X[,which(grp==g),drop=FALSE]))
        colnames(gX) = groups
        dim(gX)
        dim(ngs$X)

        dbg("getSingleSampleEnrichment:: 2")
        
        ## for large datasets pre-grouping is faster
        ss.bygroup  <- calcSingleSampleValues(gX, y, method=c("rho","gsva"))
        do.rho   = TRUE
        do.gsva  = FALSE
        do.exact = FALSE
        if(DEV.VERSION) {
            if(input$sig_ssstats=="rho") { do.rho=TRUE; do.gsva = FALSE; do.exact=TRUE }
            if(input$sig_ssstats=="gsva") { do.rho=FALSE; do.gsva = TRUE; do.exact=TRUE }
            if(input$sig_ssstats=="grp.gsva") { do.rho=FALSE; do.gsva=TRUE; do.exact=FALSE }
            if(input$sig_ssstats=="rho+gsva") { do.rho=TRUE; do.gsva=TRUE; do.exact=TRUE }
            if(input$sig_ssstats=="rho+grp.gsva") { do.rho=TRUE; do.gsva=TRUE; do.exact=FALSE }
        }
        ss.bysample <- c()
        if(do.rho) {
            ss1 <- calcSingleSampleValues(X[,], y, method=c("rho"))
            ss.bysample <- cbind(ss.bysample, rho=ss1)
        }

        dbg("getSingleSampleEnrichment:: 3")

        if(do.gsva) {
            if(do.exact) {
                ##ss.bysample <- calcSingleSampleValues(X[,1:250], y, method=c("rho","gsva"))
                ss1 <- calcSingleSampleValues(X[,], y, method=c("gsva"))
                ss.bysample <- cbind(ss.bysample, gsva=ss1)
            } else {
                ss.bysample <- calcSingleSampleValues(X, y, method=c("rho"))
                ss1 <- ss.bygroup[,"gsva"][grp]
                ss.bysample <- cbind(ss.bysample, gsva=ss1)
            }
        }

        dbg("getSingleSampleEnrichment:: done!\n")
        
        res <- list( by.sample=ss.bysample, by.group=ss.bygroup)
        return(res)
    })


    sigCalculateGSEA <- reactive({
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)
        require(fgsea)
        
        ## observe input list
        gset = head(rownames(ngs$X),100)
        gset <- getCurrentMarkers()
        if(is.null(gset)) return(NULL)
        ##if(is.null(input$sig_enplotsdb)) return(NULL)
        
        ## get all logFC of this dataset
        F <- NULL   
        ##F <- sapply(ngs$gx.meta$meta,function(x) unclass(x$fc)[,"trend.limma"])
        F <- sapply(ngs$gx.meta$meta,function(x) x$meta.fx)
        rownames(F) <- rownames(ngs$gx.meta$meta[[1]])    
        ext.db <- input$cmp_querydataset
        if(is.null(ext.db)) return(NULL)
        if(1 && ext.db!="" && ext.db!="<this dataset>") {
            if(ext.db=="<all>") {
                F <- PROFILES[["FC"]]
            } else {
                ext.db0 <- gsub("\\[|\\]","",ext.db)
                jj <- grep(ext.db0, colnames(PROFILES[["FC"]]))
                F <- PROFILES[["FC"]][,jj,drop=FALSE]
            }
        }

        ## cleanup matrix
        F = as.matrix(F)
        dim(F)
        F =  F[,which(!duplicated(colnames(F))),drop=FALSE]
        
        ## cleanup names and uppercase for mouse genes
        rownames(F) <- toupper(sub(".*:","",rownames(F)))
        gset <- toupper(sub(".*:","",gset))
        gset <- intersect(toupper(gset), rownames(F))
        length(gset)

        if(length(gset)==0) {
            cat("FATAL:: sigCalculateGSEA : gset empty!\n")
            return(NULL)
        }
        
        ## ------------ prioritize with quick correlation
        ##require(wCorr)
        y = 1*(toupper(rownames(F)) %in% toupper(gset))
        ss.rank <- function(x) scale(sign(x)*rank(abs(x)),center=FALSE)[,1]
        rho = cor(apply(F,2,ss.rank), y, use="pairwise")[,1]
        ##wt = c(mean(y==0),mean(y==1))[1+y]
        ##wt.rho = apply(F,2, function(x) weightedCorr((x), y, weights=wt, method="Pearson"))
        rho[is.na(rho)] <- 0
        names(rho) = colnames(F)
        
        ## ------- restrict to top 100 comparisons (fgsea is otherwise to
        ## ------- slow) but we do not expect many with so many
        ## ------- comparisons
        ntop = 100    
        jj <- head(order(-abs(rho)), ntop )
        F <- F[,jj,drop=FALSE]
        F <- F[!duplicated(rownames(F)),,drop=FALSE]
        F <- F + 1e-4*matrix(rnorm(length(F)),nrow(F),ncol(F))
        dim(F)
        
        ## ------------- do fast GSEA
        require(parallel)
        require(fgsea)
        gmt = list("gset"=unique(gset))
        res <- NULL
        enrich_method="rcor"
        enrich_method="fgsea"
        ##enrich_method <- input$sig_rankmethod
        
        if(enrich_method=="fgsea") {
            i=1
            withProgress(message="computing GSEA ...", value=0.8, {
                res <- lapply(1:ncol(F), function(i) {
                    res = fgsea(gmt, stats=F[,i], nperm=1000)
                    res = as.data.frame(res[,c("pval","padj","ES","NES")])
                    rownames(res)[1] = colnames(F)[i]
                    return(res)
                })
            })
            dbg("sigCalculateGSEA:: fgsea done!\n")
            res1 <- data.frame(do.call(rbind, res))
            res1$ES <- NULL
        } else {
            i=1
            fx <- 1*(rownames(F) %in% gmt[[1]])
            rho <- cor(apply(F,2,rank,na.last="keep"), fx, use="pairwise")[,1]
            pv <- cor.pvalue(rho, nrow(F))
            qv <- p.adjust(pv,method="fdr")
            res1 <- data.frame(pval=pv, padj=qv, rho=rho, NES=NA)
            rownames(res1) <- names(pv)
        }

        ## columns are: NES, pval, fdr, contrast
        res1 <- as.matrix(res1)
        res1 <- res1[match(colnames(F),rownames(res1)),,drop=FALSE]
        
        if( nrow(res1) != ncol(F)) {
            cat("WARNING sigCalculateGSEA:: fgsea results are corrupted?\n")
            cat("WARNING sigCalculateGSEA:: got contrasts: ",res$contrast,"\n")
            cat("WARNING sigCalculateGSEA:: colnames.F= ",colnames(F),"\n")
        }
        
        ## make nice table
        ##nes   <- unlist(sapply(res, function(x) x$NES))
        ##pval  <- unlist(sapply(res, function(x) x$pval))
        nes <- res1[,"NES"]
        pval <- res1[,"pval"]
        qval <- p.adjust( pval, method="fdr")
        rho <- rho[colnames(F)]
        
        output <- as.matrix(cbind(NES=nes, p=pval, q=qval, rho=rho))
        rownames(output) <- colnames(F)    
        output <- output[order(-abs(output[,"NES"])),,drop=FALSE]
        F <- F[,rownames(output),drop=FALSE]    
        gsea <- list(F=as.matrix(F), gset=gset, output=output)
        dbg("sigCalculateGSEA:: done!\n")
        return(gsea)
    })

    ##X=gX;method=c("rho","gsva") 
    calcSingleSampleValues <- function(X, y, method=c("rho","gsva") ) {
        ##
        ## Calculates single-sample enrichment values for given matrix and
        ## binarized signature vector.
        ##
        ##
        ## very fast rank difference

        dbg("<signature:calcSingleSampleValues> called\n")
        
        if(is.null(names(y)) && length(y)!=nrow(X) ) {
            cat("<signature:calcSingleSampleValues> FATAL ERROR: y must be named if not matched\n")
            return(NULL)
        }
        
        if(!is.null(names(y)) && length(y)!=nrow(X) ) {
            y <- y[match(rownames(X),names(y))]
        }
        names(y) <- rownames(X)
        jj <- which(!is.na(y))
        X <- X[jj,]
        y <- y[jj]

        dbg("<signature:calcSingleSampleValues> 1\n")
        
        if(sum(y!=0)==0) {
            cat("<signature:calcSingleSampleValues> WARNING: y is all zero!\n")        
            matzero <- matrix(0, nrow=ncol(X), ncol=length(method))
            colnames(matzero) <- method
            rownames(matzero) <- colnames(X)
            return(matzero)
        }
        ss.rank <- function(x) scale(sign(x)*rank(abs(x)),center=FALSE)[,1]
        
        S = list()
        if("rho" %in% method) {
            S[["rho"]] <- cor(apply(X, 2, ss.rank), y, use="pairwise")[,1]
            ##S$rho <- scale(S$rho)[,1]  ## should we scale??
        }

        dbg("<signature:calcSingleSampleValues> 2\n")
        
        ## calculate GSVA
        if("gsva" %in% method) {
            require(GSVA)
            require(parallel)
            gset = names(y)[which(y!=0)]
            gmt <- list("gmt"=gset)
            res.gsva <- GSVA::gsva( X, gmt, method="gsva", parallel.sz=1) ## parallel=buggy
            res.colnames = colnames(res.gsva)
            fc = as.vector(res.gsva[1,])
            names(fc) = res.colnames
            S[["gsva"]] = fc[colnames(X)]
        }    
        s.names = names(S)
        if(length(S)>1) {
            S1 = do.call(cbind, S)
        } else {
            S1 <- S[[1]]
        }

        dbg("<signature:calcSingleSampleValues> done!\n")

        S1 = as.matrix(S1)
        rownames(S1) = colnames(X)
        colnames(S1) = s.names
        ##S1 <- S1[order(-S1[,"gsva"]),]
        return(S1)
    }


    ##================================================================================
    ## Enrichment {data-height=800}
    ##================================================================================
    
    sig_enplotsFUNC <- reactive({
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)
        
        ##require(shinycssloaders)
        gsea <- sigCalculateGSEA()
        if(is.null(gsea)) return(NULL)
        
        F <- as.matrix(gsea$F)
        gset <- gsea$gset
        qv <- gsea$output[,"q"]
        
        require(gplots)
        cex.main=1.1
        par(mfrow=c(5,3), mar=c(0.3,4,3,1), mgp=c(2.2,0.8,0) )
        if(ncol(F)>15) {
            par(mfrow=c(6,5), mar=c(0.2,2,3,0.6))
            cex.main=0.9
        }
        ## if(ncol(F)>24) par(mfrow=c(7,5), mar=c(1,2,2.5,0.6))
        for(i in 1:min(30,ncol(F))) {
            f <- colnames(F)[i]
            tt <- sub(".*\\]","",f)
            tt <- breakstring(substring(tt,1,50),28,force=TRUE)
            gsea.enplot(F[,i], gset, main=tt, cex.main=cex.main)
            qv1 <- paste("q=",round(qv[i],digits=4))
            legend("topright",qv1, cex=0.85, bty="n", adj=0)
            if(grepl("^\\[",f)) {
                db <- sub("\\].*","]",colnames(F)[i])
                legend("topleft",db, cex=0.85, bty="n", adj=0)
            }
        }
    })

    observe({
        cmapsets <- c(sort(unique(gsub("\\].*","]",colnames(PROFILES$FC)))))
        cmapsets <- c("<this dataset>",cmapsets,"<all>")
        ##updateSelectInput(session, "cmp_cmapsets1", choices=cmapsets, selected="<this>")
        updateSelectInput(session, "cmp_querydataset", choices=cmapsets)
    })


    sig_enplots_info = "The <strong>Enrichment</strong> tab performs the enrichment analysis of the gene list against all contrasts by running the GSEA algorithm and plots enrichment outputs. Enrichment statistics can be found in the corresponding table"
    
    sig_enplots_caption = "<b>Enrichment plots.</b> This figure shows the enrichment of the query signature in all constrasts. Positive enrichment means that this particular contrast shows similar expression changes as the query signature."

    sig_enplots.opts = NULL
    sig_enplots.module <- plotModule(
        "sig_enplots", sig_enplotsFUNC, plotlib="base",
        info.text = sig_enplots_info,
        options = sig_enplots.opts,
        pdf.width=8, pdf.height=8, res=90
    )
    output <- attachModule(output, sig_enplots.module)

    output$sig_enplots_UI <- renderUI({
        fillCol(
            height = 750,
            flex = c(1, NA),
            moduleWidget(sig_enplots.module, ns=ns),
            div(HTML(sig_enplots_caption), class="caption")
        )
    })
    

    ##================================================================================
    ## Overlap/similarity
    ##================================================================================

    sig_getOverlapTable <- reactive({
        ##
        ##
        ##
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)
        
        markers <- GSETS[[100]]
        markers <- getCurrentMarkers()
        if(is.null(markers)) return(NULL)
        
        ## fold change just for ranking of genes
        F <- sapply(ngs$gx.meta$meta, function(x) unclass(x$fc)[,"trend.limma"])
        fx <- rowMeans(F**2)
        ##markers=head(names(sort(-abs(fx))),100)
        
        ## fisher test
        ##ii <- setdiff(match(markers, colnames(GSETxGENE)),NA)
        ii <- setdiff(match(toupper(markers), colnames(GSETxGENE)),NA)
        N <- cbind(k1=Matrix::rowSums(GSETxGENE!=0), n1=ncol(GSETxGENE),
                   k2=Matrix::rowSums(GSETxGENE[,ii]!=0), n2=length(ii) )
        rownames(N) = rownames(GSETxGENE)
        ##N <- N[which(!(N[,1]==0 & N[,3]==0)), ]
        N <- N[which(N[,1]>0 | N[,3]>0), ]
        odds.ratio = ( N[,3]/ N[,4]) / ( N[,1]/ N[,2]) 
        dim(N)
        
        require(corpora) 
        dbg("computing Fisher-test p values...")
        pv <- corpora::fisher.pval( N[,1], N[,2], N[,3], N[,4], log.p=FALSE)
        dbg("done!\n")
        head(pv)
        names(pv) <- rownames(N)
        pv = pv[match(names(odds.ratio),names(pv))]
        ##qv = p.adjust(pv, method="fdr")
        qv = p.adjust(pv, method="bonferroni")
        A = data.frame( odds.ratio=odds.ratio, p.fisher=pv, q.fisher=qv)
        dim(A)

        ## limit the list??
        table( qv < 0.05)
        table( qv < 0.2)
        table( qv < 0.999)
        A <- A[which( A$q.fisher < 0.999),]
        ##A <- A[which( A$q.fisher < 0.05),]
        dim(A)
        
        if(0) {
            ## Gene set clustering so that similar gene sets are shown
            ## together??? Bit like Gprofiler
            ##
            require(qlcMatrix)
            require(nclust)
            G = GSETxGENE[rownames(A),]
            G = G[which(Matrix::rowSums(G!=0) > 10), ] ## at least 10 genes
###G = G[, which(Matrix::colSums(G!=0) >= 10) ] ## at least 10 genesets
            dim(G)        
            hc <- fastcluster::hclust(as.dist(cosSparse(t(G))))
            ##hc <- nclust(as.matrix(G))
            A <- A[rownames(G)[hc$order],]
        }

        ## get shared genes
        dbg("determining shared genes...\n")
        aa = rownames(A)
        y <- 1*(colnames(GSETxGENE) %in% toupper(markers))
        names(y) <- colnames(GSETxGENE)
        ncommon <- Matrix::colSums(t(GSETxGENE[aa,])*as.vector(y)!=0)
        ntotal  <- Matrix::rowSums(GSETxGENE[aa,]!=0)
        A$ratio <- ncommon / ntotal
        ratio.kk <- paste0(ncommon,"/",ntotal)    
        
        gg <- colnames(GSETxGENE)
        gset <- names(y)[which(y!=0)]
        G1 = GSETxGENE[aa,which(y!=0)]
        commongenes <- apply(G1, 1, function(x) colnames(G1)[which(x!=0)])
        ##commongenes <- lapply(commongenes, function(x) x[order(-fx[x])])
        ##commongenes <- mclapply(commongenes, function(x) x[order(-fx[x])])
        for(i in 1:length(commongenes)) {
            gg <- commongenes[[i]]
            gg <- gg[order(-abs(fx[gg]))]
            if(length(gg)>10) {
                others <- paste0("(+",length(gg)-10," others)")
                gg <- c(head(gg,10),others)
            }
            commongenes[[i]] <- paste(gg,collapse=",")        
        }
        ##commongenes <- sapply(commongenes,paste,collapse=",")        
        commongenes <- unlist(commongenes)

        ## construct results dataframe
        gset.names <- substring(rownames(A),1,72)    
        ##aa <- apply(A, 2, formatC, format="e", digits=3)
        A$ratio <- round(A$ratio, digits=3)
        A$log.OR <- round(log10(A$odds.ratio), digits=3)
        A$odds.ratio <- round(A$odds.ratio, digits=3)
        db = sub(":.*","",gset.names)
        score = (log10(A$odds.ratio) * -log10(A$q.fisher + 1e-40))**0.5
        score = round(score, digits=3)
        df <- cbind(db=db, geneset=gset.names, score=score, "k/K"=ratio.kk, A, common.genes=commongenes)
        if(DEV.VERSION) {
            df <- df[,c("db","geneset","score","k/K","ratio","odds.ratio","log.OR","q.fisher","common.genes")]
        } else {
            df <- df[,c("db","geneset","score","k/K","odds.ratio","q.fisher","common.genes")]
        }
        
        ##df <- df[order(-df$odds.ratio),]
        df <- df[order(-df$score),]
        return(df)
    })

    sig_overlapScorePlot.RENDER <- reactive({
        
        df <- sig_getOverlapTable()
        req(df)
        req(input$sig_overlapTable_rows_all)

        sel <- 1:nrow(df)
        sel<- input$sig_overlapTable_rows_all
        df1 <- df[sel,]
        df1$geneset = as.character(rownames(df1))
        df1$db = factor(df1$db)
        
        ntop = 1000
        ntop = as.integer(input$sig_overlapScorePlot_ntop)
        df1 = df1[head(order(-df1$score),ntop),]
        jj = order(df1$db, -df1$score)
        df1 = df1[jj,]

        df1$idx = factor(1:nrow(df1), levels=1:nrow(df1))
        klr = rep(brewer.pal(8,"Set2"),10)[as.integer(df1$db)]
        
        plt <- plot_ly(
            df1, x = ~idx, y = ~score,
            type='bar',  ## orientation='v',
            text = ~geneset, hoverinfo = 'text',
            hovertemplate = paste0("%{text}<br>score: %{y}<extra>",df1$db,"</extra>"),
            ##hovertemplate = "%{y}",
            marker = list( color=klr ) ) %>%
            layout(
                showlegend = FALSE,
                dragmode= 'select',
                ##annotations = anntitle(colnames(rho)[i]),
                ##annotations = list(text="TITLE"),
                yaxis = list(##range = c(0,1),
                    titlefont = list(size=11),
                    tickfont = list(size=10),
                    showgrid = TRUE,
                    title = "overlap score" ),
                xaxis = list(
                    title = "",
                    showgrid = FALSE,
                    showline = FALSE,
                    showticklabels = FALSE,
                    showgrid = FALSE,
                    zeroline = FALSE)) 

        if( min(nrow(df1),ntop) < 100 && input$sig_overlapScorePlot_shownames) {
            ## labeling the y-axis inside bars
            plt <- plt %>%
                add_annotations( yref='paper', xref = 'x',
                                x = ~idx, y=0.005, yanchor='bottom',
                                text = substring(df1$geneset,1,35),
                                textangle = -90,
                                font = list(size = 10),
                                showarrow = FALSE, align='right')
        }

        plt    
    })

    sig_overlapTable.RENDER <- reactive({

        df <- sig_getOverlapTable()
        req(df)    

        df$geneset <- wrapHyperLink(df$geneset, df$geneset)

        numeric.cols <- which(sapply(df, is.numeric))
        numeric.cols <- intersect(c("p.fisher","q.fisher"),colnames(df))
        
        DT::datatable(df, class='compact cell-border stripe',
                      rownames=FALSE, escape = c(-1,-2),
                      extensions = c('Scroller'),
                      selection='none',
                      fillContainer=TRUE,
                      options=list(
                          dom = 'lfrtip',
                          ## pageLength = 40, ##lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE, scrollY = 500, scroller=TRUE ## deferRender=TRUE,
                      )  ## end of options.list 
                      ) %>%
            formatSignif(numeric.cols,4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
                DT::formatStyle("score",
                                background = color_from_middle( df$score, 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')

    })

    sig_overlapScorePlot.opts = tagList(
        tipify(radioButtons(ns("sig_overlapScorePlot_ntop"),"Number of features",c(60,120,250),inline=TRUE),
               "Specify the number to top features to show.", placement="top", options = list(container = "body")),
        tipify(checkboxInput(ns("sig_overlapScorePlot_shownames"),"Show feature names",TRUE),
               "Select to show/hide the feature names in the plot.", placement="top", options = list(container = "body"))
    )


    sig_overlapScorePlot.module <- plotModule(
        "sig_overlapScorePlot", sig_overlapScorePlot.RENDER, plotlib="plotly",
        title = "Signature overlap scores", label="a",
        info.text = "Top overlapping gene sets with selected signature. The vertical axis shows the overlap score of the gene set which combines the odds ratio and significance (q-value) of the Fisher's test.",
        options = sig_overlapScorePlot.opts,
        pdf.width=12, pdf.height=6, res=100
    )
    output <- attachModule(output, sig_overlapScorePlot.module)

    
    sig_overlapTable.module <- tableModule(
        "sig_overlapTable", sig_overlapTable.RENDER,
        title = "Overlap with other signatures", label="b",
        ## just.info=TRUE, no.download=TRUE,
        info.text = "Under the <strong>Overlap/similarity tab</strong>, users can find the similarity of their gene list with all the gene sets and pathways in the platform, including statistics such as the total number of genes in the gene set (K), the number of intersecting genes between the list and the gene set (k), the overlapping ratio of k/K, logarithm of the  odds ratio (log.OR), as well as the p and q values by the Fisher’s test for the overlap test.",
        ##options = sig_overlapTable.opts,
        )
    output <- attachModule(output, sig_overlapTable.module)


    sig_overlap_caption = "<b>Overlap/Similarity table.</b><b>(a)</b> Top overlapping gene sets with selected signature. The vertical axis shows the overlap score of the gene set which combines the odds ratio and significance (q-value) of the Fisher's test. <b>(b)</b> Table summarizing the results of the Fishers's test for overlap. The column \'common genes\' reports the shared gene in order of largest fold-change."
    

    output$sig_overlapAnalysis_UI <- renderUI({
        fillCol(
            flex = c(1,0.04,1,0.08,NA),
            height = 750,
            moduleWidget(sig_overlapScorePlot.module, outputFunc="plotlyOutput", ns=ns),
            br(),
            moduleWidget(sig_overlapTable.module, outputFunc="dataTableOutput", ns=ns),
            br(),
            div(HTML(sig_overlap_caption), class="caption")
        )
    })


    ##================================================================================
    ## Markers {data-height=800}
    ##================================================================================

    sig_markersplotFUNC <- reactive({
        ##if(!input$tsne.all) return(NULL)
        require(RColorBrewer)
        require(gplots)
        
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)

        dbg("<signature:sig_markersplotFUNC> called\n")        
        
        markers <- ngs$families[[2]]
        markers <- COLLECTIONS[[10]]
        markers <- getCurrentMarkers()
        if(is.null(markers)) return(NULL)
        
        level = "gene"
        ##markers <- intersect(markers,ngs$genes$gene_name)
        ##jj <- match(markers,ngs$genes$gene_name)
        xgene <- ngs$genes[rownames(ngs$X),]$gene_name
        jj <- match(toupper(markers), toupper(xgene))
        jj <- setdiff(jj,NA)
        gx <- ngs$X[jj,,drop=FALSE]
        
        if(nrow(gx)==0) {
            cat("WARNING:: Markers:: markers do not match!!\n")
            return(NULL)
        }
        
        ## get t-SNE positions of samples
        pos = ngs$tsne2d[colnames(gx),]
        gx = gx - min(gx,na.rm=TRUE) + 0.001 ## subtract background
        dim(gx)
        grp <- ngs$samples[colnames(gx),"group"]
        zx <- t(apply(gx,1,function(x) tapply(x,as.character(grp),mean)))
        gx <- gx[order(-apply(zx,1,sd)),,drop=FALSE]
        rownames(gx) = sub(".*:","",rownames(gx))
        
        ## ---------------- get GSVA values
        res <- getSingleSampleEnrichment()
        if(is.null(res)) return(NULL)
        
        S <- res$by.sample
        if(NCOL(S)==1) {
            fc = S[,1]
        } else {
            fc = colMeans( t(S) / (1e-8+sqrt(colSums(S**2))) ) ## scaled mean
        }
        fc <- scale(fc)[,1]  ## scale??
        names(fc) = rownames(S)
        ##fc1 = sign(fc) * (fc/(1e-8+max(abs(fc))))**2
        fc1 = tanh(1.0*fc / (1e-4+sd(fc)))
        fc1 = fc1[rownames(pos)]
        
        cex1 = 1.1
        if(nrow(pos) < 30) cex1 = 1.8
        if(nrow(pos) > 200) cex1 = 0.5

        cex1 = 1.2
        cex1 <- 0.85*c(1.6,1.2,0.8,0.5)[cut(nrow(pos),breaks=c(-1,40,200,1000,1e10))]    
        cex2 <- ifelse(level=="gene",1,0.8)
        klrpal = colorRampPalette(c("grey90", "grey60", "red3"))(16)

        nmax = 35
        par(mfrow=c(6,6), mar=c(0,0.2,0.5,0.2), oma=c(2,1,2,1)*0.8 )
        top.gx = head(gx,nmax)
        top.gx = top.gx[order(rownames(top.gx)),,drop=FALSE]
        
        i=1    
        for(i in 0:min(nmax,nrow(top.gx))) {
            jj <- 1:ncol(top.gx)
            if(i==0) {
                klr1 = BLUERED(16)[8 + round(7*fc1)]
                tt = "INPUT SIGNATURE"
                jj <- order(abs(fc1))
            } else {
                colvar = pmax(top.gx[i,],0) 
                colvar = 1+round(15*(colvar/(0.7*max(colvar)+0.3*max(top.gx))))
                klr1 = klrpal[colvar]
                gene <- substring(sub(".*:","",rownames(top.gx)[i]),1,80)
                tt <- breakstring(gene, n=20, force=TRUE)
                jj <- order(abs(top.gx[i,]))
            }
            klr1 = paste0(col2hex(klr1),"99")
            
            plot( pos[jj,], pch=19, cex=cex1, col=klr1[jj],
                 xlim=1.2*range(pos[,1]), ylim=1.2*range(pos[,2]),
                 fg = gray(ifelse(i==0,0.1,0.8)), bty = "o",
                 xaxt='n', yaxt='n', xlab="tSNE1", ylab="tSNE2")
            legend("topleft", tt, cex=cex2, col="grey30", text.font=ifelse(i==0,2,1),
                   inset=c(-0.1,-0.05), bty="n")

        }

        dbg("<signature:sig_markersplotFUNC> done!\n")        
        
    })

    ##sig_markers.opts = tagList()
    sig_markers.module <- plotModule(
        "sig_markersplot", sig_markersplotFUNC, plotlib="base",
        info.text = "After uploading a gene list, the <strong>Markers</strong> section produces a t-SNE plot of samples for each gene, where the samples are colored with respect to the upregulation (in red) or downregulation (in blue) of that particular gene.",
        ##options = sig_markers.opts,
        pdf.width=8, pdf.height=8, res=100
    )
    output <- attachModule(output, sig_markers.module)

    sig_markers_caption = "<b>Markers t-SNE plot</b>. This figure shows the t-SNE plot for each gene, where the dot (corresponding to samples) are colored depending on the upregulation (in red) or downregulation (in blue) of that particular gene."

    output$sig_markers_UI <- renderUI({
        fillCol(
            flex = c(1,0.04,NA),
            height = 750,
            moduleWidget(sig_markers.module, ns=ns),
            br(),
            div(HTML(sig_overlap_caption), class="caption")
        )
    })
    
    ##================================================================================
    ## Enrichment {data-height=800}
    ##================================================================================
    
    sig_enrichmentByContrastTable.RENDER <- reactive({
        
        gsea <- sigCalculateGSEA()
        if(is.null(gsea)) return(NULL)
        
        output <- as.matrix(gsea$output)
        output <- round(output, digits=4)
        output <- data.frame( contrast=rownames(output), output)
        if(!DEV.VERSION) {
            output$p <- NULL
            output$rho <- NULL
        }
        
        color_fx = as.numeric(output[,"NES"])
        color_fx[is.na(color_fx)] <- 0  ## yikes...    
        numeric.cols <- which(sapply(output, is.numeric))
        numeric.cols

        DT::datatable(output, class='compact cell-border stripe',
                      rownames=FALSE,
                      ##extensions = c('Buttons','Scroller'),
                      extensions = c('Scroller'),
                      ##selection='none',
                      selection=list(mode='single', target='row', selected=c(1)),
                      options=list(
                          ##dom = 'Blfrtip',
                          dom = 'lrtip',
                          ## pageLength = 40, ##lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE, scrollY = 220, scroller=TRUE, deferRender=FALSE,
                          buttons = c('copy','csv','pdf')
                      )) %>%  ## end of options.list 
            formatSignif(numeric.cols,4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
                DT::formatStyle("NES",
                                background = color_from_middle(color_fx, 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')

    })

    sig_enrichmentByContrastGenes.RENDER <- reactive({
        
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)

        gsea <- sigCalculateGSEA()
        if(is.null(gsea)) return(NULL)
        
        i=1
        i <- input$sig_enrichmentByContrastTable_rows_selected
        if(is.null(i) || length(i)==0) return(NULL)
        contr <- rownames(gsea$output)[i]

        ext.db <- input$cmp_querydataset
        if(is.null(ext.db)) return(NULL)

        ##if(input$sig_enplotsdb=="this dataset") {
        if(ext.db=="<this dataset>") {
            fc <- ngs$gx.meta$meta[[contr]]$meta.fx
            qv <- ngs$gx.meta$meta[[contr]]$meta.q
            names(fc) <- rownames(ngs$gx.meta$meta[[contr]])
            names(qv) <- rownames(ngs$gx.meta$meta[[contr]])
            names(fc) <- toupper(names(fc))
            names(qv) <- toupper(names(qv))
        } else {
            ##load(file.path(FILES,"allFoldChanges-pub.rda"))
            fc <- PROFILES[["FC"]][,contr]
            qv <- rep(NA,length(fc))
            names(qv) <- names(fc)
        }

        gset <- getCurrentMarkers()
        if(is.null(gset)) return(NULL)
        
        gset <- setdiff(toupper(gset),c("",NA))
        genes <- intersect(gset,names(fc))    
        fc <- fc[genes]
        fc <- round(fc[order(-abs(fc))],digits=3)
        qv <- qv[names(fc)]
        gene.tt <- substring(GENE.TITLE[toupper(names(fc))],1,40)
        names(gene.tt) <- names(qv) <- names(fc)    
        output <- data.frame(gene=names(fc), title=gene.tt, FC=fc, q=qv)
        
        color_fx = as.numeric(output[,"FC"])
        color_fx[is.na(color_fx)] <- 0  ## yikes...

        numeric.cols <- which(sapply(output, is.numeric))
        numeric.cols
        
        DT::datatable(output, class='compact cell-border stripe',
                      rownames=FALSE,
                      ##extensions = c('Buttons','Scroller'),
                      extensions = c('Scroller'), selection='none',
                      options=list(
                          ##dom = 'Blfrtip',
                          dom = 'lrftip',
                          ## pageLength = 40, ##lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE, scrollY = 210, scroller=TRUE, deferRender=FALSE,
                          buttons = c('copy','csv','pdf')
                      )) %>%  ## end of options.list 
            formatSignif(numeric.cols,4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
                DT::formatStyle("FC",
                                background = color_from_middle(color_fx, 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')    
    })


    info.text1 = "The table summarizes the enrichment statistics of the gene list against all contrasts by running the GSEA algorithm and plots enrichment outputs."
    sig_enrichmentByContrastTable_module <- tableModule(
        id = "sig_enrichmentByContrastTable", func = sig_enrichmentByContrastTable.RENDER,
        info.text = info.text1,
        title = "Enrichment by contrasts", label="a"
    )
    output <- attachModule(output, sig_enrichmentByContrastTable_module)


    info.text2 = "This table shows the genes of the current signature."
    sig_enrichmentByContrastGenes_module <- tableModule(
        id = "sig_enrichmentByContrastGenes", func = sig_enrichmentByContrastGenes.RENDER,
        info.text = info.text2,
        title = "Genes in signature", label="b"
    )
    output <- attachModule(output, sig_enrichmentByContrastGenes_module)
    
    sig_enrichmentTables_caption = "<b>Enrichment of query signature across all contrasts.</b> <b>(a)</b> Enrichment scores across all contrasts for the selected query signature . The NES corresponds to the normalized enrichment score of the GSEA analysis.  <b>(b)</b> Genes in the query signature sorted by decreasing (absolute) fold-change corresponding to the contrast selected in Table (a)."

    output$sig_enrichmentTables_UI <- renderUI({
        fillCol(
            flex = c(1.0,0.04,1.1,0.04,NA), ## width = 600,
            height = 750,
            moduleWidget(sig_enrichmentByContrastTable_module, outputFunc="dataTableOutput"),
            br(),
            moduleWidget(sig_enrichmentByContrastGenes_module, outputFunc="dataTableOutput"),
            br(),
            div(HTML(sig_enrichmentTables_caption), class="caption")
        )
    })

}
