##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##=============================================================================
##==========    Platform global and user settings =============================
##=============================================================================

USER.GENESETTEST.METHODS <- NULL
USER.GENETEST.METHODS <- NULL

##=============================================================================
##===================    Platform helper functions ============================
##=============================================================================

pgx.phenoMatrix <- function(pgx, phenotype) {
    y <- pgx$samples[,phenotype]
    mm <- t(model.matrix(~0+y))
    rownames(mm) <- sub("^y","",rownames(mm))
    colnames(mm) <- rownames(pgx$samples)
    as.matrix(mm)
}

cex=1
text_repel.NOTWORKING <- function( x, y, text, cex=1, force=1e-7, maxiter=20000)
{    
    if(0) {
        ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg, label=rownames(mtcars))) +
            ggplot2::geom_point() +
            ggrepel::geom_text_repel()
        
        x=mtcars[,"wt"]
        y=mtcars[,"mpg"]
        text=rownames(mtcars)
        
        labx <- out[,3]
        laby <- out[,4]
        ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
            ggplot2::geom_point() +
            ggplot2::geom_text(x = labx, y = laby, label=rownames(mtcars))
                      

    }
    ## x and y posiitons as a dataframe
    df <- data.frame(x=x, y=y, text=text)
    w <- diff(range(x))
    h <- diff(range(y))
    dx0 <- w * 0.08
    dy0 <- h * 0.05
    par(mfrow=c(1,1)); plot(x,y,type='n')
    dx1 <- max(strwidth(text, cex=cex, units="user") * w)
    dy1 <- max(strheight(text, cex=cex, units="user") * w)
    out <- util.findboxes(
        df, "x", "y",
        box_padding_x = dx1,
        box_padding_y = dy1,
        point_padding_x = dx0,
        point_padding_y = dy0,
        xlim = range(x),
        ylim = range(y),
        force = 1e-7, maxiter = 20000
    )
    out[,3:4]
}

pos.compact <- function(pos, d=0.01) {
    ## make positions more dense removing white space
    for(i in 1:ncol(pos)) {
        x=pos[,i]
        dr = d*diff(range(x))
        names(x)=1:nrow(pos)
        ii <- order(x)
        x1 = cumsum(c(x[ii[1]],pmin(diff(x[ii]),dr)))
        pos[,i] = x1[order(as.integer(names(x1)))]
    }
    pos
}

#' Given a Set of Points and Box sizes,
#' https://github.com/slowkow/ggrepel/issues/24
util.findboxes <- function( df, xcol, ycol,
                           box_padding_x, box_padding_y,
                           point_padding_x, point_padding_y,
                           xlim, ylim,
                           force = 1e-7, maxiter = 20000
                           )
{
    
    
  # x and y posiitons as a dataframe
  posdf <- df[c(xcol, ycol)]

  # returnd a df where columns are points
  boxdf <- apply(posdf, 1, function(row) {
    xval <- row[xcol]
    yval <- row[ycol]
    return(c(
      xval - box_padding_x / 2,
      yval - box_padding_y / 2,
      xval + box_padding_x / 2,
      yval + box_padding_y / 2
    ))
  })
  # columns are x1,y1,x2,y2
  boxmatrix <- as.matrix(t(boxdf))

  moved <- ggrepel:::repel_boxes(
    data_points = as.matrix(posdf),
    point_padding_x = point_padding_x,
    point_padding_y = point_padding_y,
    boxes = boxmatrix,
    xlim = xlim,
    ylim = ylim,
    hjust = 0.5,
    vjust = 0.5,
    force = force,
    maxiter = maxiter
  )

  finaldf <- cbind(posdf, moved)
  names(finaldf) <- c("x1", "y1", "x2", "y2")
  return(finaldf)
}

star.symbols <- function(n, pch="\u2605") {
    if(n==0) return("")
    paste(rep(pch,n),collapse="")
}

search_path <- function(paths, file) {
    dir <- paths[which(file.exists(file.path(paths,file)))]
    if(length(dir)==0) return(NULL)
    file.path(dir[1],file)
}

rowscale <- function(x) {
    x  <- x - Matrix::rowMeans(x,na.rm=TRUE)
    x / (1e-4 + sqrt(rowMeans(x**2,na.rm=TRUE)))
}

strwrap2 <- function(str,n ) {
    sapply(str,function(s) paste(base::strwrap(s,n),collapse="\n"))
}

add_opacity <- function(hexcol,opacity) {
    ##toRGB(hexcol)
    col1 <- rep(NA,length(hexcol))
    ii <- which(!is.na(hexcol))
    rgba <- strsplit(gsub("rgba\\(|\\)","",plotly::toRGB(hexcol[ii],opacity)),split=",")   
    rgba <- apply(do.call(rbind, rgba),2,as.numeric)
    if(length(hexcol)==1) rgba <- matrix(rgba,nrow=1)
    col1[ii] <- rgb(rgba[,1]/255,rgba[,2]/255,rgba[,3]/255,rgba[,4])
    col1
}

logCPM <- function(counts, total=1e6, prior=1) {
    ## Transform to logCPM (log count-per-million) if total counts is
    ## larger than 1e6, otherwise scale to previous avarage total count.
    ##
    ##
    if(is.null(total)) {
        ##total <- nrow(counts)
        ##total <- mean(colSums(counts1!=0)) ## avg. number of expr genes
        total0 <- mean(Matrix::colSums(counts,na.rm=TRUE)) ## previous sum 
        total <- ifelse( total0 < 1e6, total0, 1e6 )
        message("[logCPM] setting column sums to = ",round(total,2))
    }
    if(any(class(counts)=="dgCMatrix")) {
        ## fast/sparse calculate CPM
        cpm <- counts
        cpm[is.na(cpm)] <- 0  ## OK??
        cpm@x <- total * cpm@x / rep.int(Matrix::colSums(cpm), diff(cpm@p))  ## fast divide by columns sum
        cpm@x <- log2(prior + cpm@x)
        return(cpm)
    } else {
        cpm <- t(t(counts) / Matrix::colSums(counts,na.rm=TRUE)) * total
        x <- log2(prior + cpm)
        return(x)
    }
}

pgx.checkObject <- function(pgx) {
    must.have <- c("counts","samples","genes","model.parameters",
                   "X","gx.meta","gset.meta","gsetX","GMT")
    not.present <- setdiff(must.have,names(pgx))
    if(length(not.present)>0) {
        not.present <- paste(not.present, collapse=" ")
        message("[pgx.checkObject] WARNING!!! object does not have: ",not.present)
    }
    all(must.have %in% names(pgx))
}

matGroupMeans <- function(X, group, FUN=rowMeans, dir=1) {
    if(dir==2) X <- t(X)
    mX <- do.call(cbind, tapply(1:ncol(X),group,function(i) rowMeans(X[,i,drop=FALSE],na.rm=TRUE)))
    if(dir==2) mX <- t(mX)
    mX
}

knnMedianFilter <- function(x, pos, k=10)
{
    
    nb <- FNN::get.knn(pos[,], k=k)$nn.index
    fx <- factor(x)
    mx <- matrix(fx[as.vector(nb)],nrow=nrow(nb),ncol=ncol(nb))
    x1 <- apply(mx,1,function(x) names(which.max(table(x))))
    x1
}

nmfImpute <- function(x,k=5) {
    ## Impute missing values with NMF
    ##
    
    k = min(k, dim(x))
    nmf <- NNLM::nnmf(x, k=k, check.k=FALSE, rel.tol = 1e-2, verbose=0)
    xhat <-  with(nmf, W %*% H);
    x[is.na(x)] <- xhat[is.na(x)]
    if(sum(is.na(x))>0) {
        nmf1 <- NNLM::nnmf(x, k=1, check.k=FALSE, rel.tol = 1e-2, verbose=0)    
        xhat1 <-  with(nmf1, W %*% H);    
        x[is.na(x)] <- xhat1[is.na(x)]
    }
    x
}

knnImputeMissing <- function(x, pos, missing=NA, k=10)
{
    
    k0 <- which(x==missing)    
    k1 <- which(x!=missing)
    if(length(k0)==0) {
        return(x)
    }
    pos0 <- pos[k0,]
    pos1 <- pos[k1,]
    nb <- FNN::get.knnx(pos1, pos0, k=k)$nn.index
    fx <- factor(x[k1])
    mx <- matrix(fx[as.vector(nb)],nrow=nrow(nb),ncol=ncol(nb))
    x.imp <- apply(mx,1,function(x) names(which.max(table(x))))
    x[which(x==missing)] <- x.imp
    x
}

randomImputeMissing <- function(x) {
    i=1
    for(i in 1:ncol(x)) {
        jj <- which(is.na(x[,i]) | x[,i]=="NA")
        if(length(jj)) {
            rr <- sample( x[-jj,i], length(jj), replace=TRUE)
            x[jj,i] <- rr
        }
    }
    return(x)
}

human2mouse.SLLOWWW <- function(x){
    
    human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genesV2 = biomaRt::getLDS(attributes = c("hgnc_symbol"),
                     filters = "hgnc_symbol",
                     values = x , mart = human,
                     attributesL = c("mgi_symbol"),
                     martL = mouse,
                     uniqueRows=T)
    genesx <- unique(genesV2[, 2])
    ## Print the first 6 genes found to the screen
    print(Matrix::head(genesx))
    return(genesx)
}

human2mouse <- function(x) {
    
    homologene::human2mouse(x)
}
mouse2human <- function(x) {
    
    homologene::mouse2human(x)
}

##type=NULL;org="human";keep.na=FALSE
probe2symbol <- function(probes, type=NULL, org="human", keep.na=FALSE)
{   
    require(org.Hs.eg.db)
    require(org.Mm.eg.db)

    ## strip postfix for ensemble codes
    if(mean(grepl("^ENS",probes))>0.5) {
        probes <- gsub("[.].*","",probes)
    }
    
    if(is.null(type)) {
        
        hs.list <- list(
            "human.ensembl" = unlist(as.list(org.Hs.egENSEMBL)),
            "human.ensemblTRANS" = unlist(as.list(org.Hs.egENSEMBLTRANS)),            
            #"human.unigene" = unlist(as.list(org.Hs.egUNIGENE)),
            "human.refseq"  = unlist(as.list(org.Hs.egREFSEQ)),
            "human.accnum"  = unlist(as.list(org.Hs.egACCNUM)),
            "human.uniprot" = unlist(as.list(org.Hs.egUNIPROT)),
            "human.symbol"  = unlist(as.list(org.Hs.egSYMBOL))
            )
        
        mm.list <- list(
            "mouse.ensembl" = unlist(as.list(org.Mm.egENSEMBL)),
            "mouse.ensemblTRANS" = unlist(as.list(org.Mm.egENSEMBLTRANS)),            
            #"mouse.unigene" = unlist(as.list(org.Mm.egUNIGENE)),
            "mouse.refseq"  = unlist(as.list(org.Mm.egREFSEQ)),
            "mouse.accnum"  = unlist(as.list(org.Mm.egACCNUM)),
            "mouse.uniprot" = unlist(as.list(org.Mm.egUNIPROT)),
            "mouse.symbol"  = unlist(as.list(org.Mm.egSYMBOL))
        )
        id.list <- c(hs.list, mm.list)
        mx <- sapply(id.list, function(id) mean(probes %in% id))
        mx
        org=type=NULL
        max.mx <- max(mx,na.rm=TRUE)
        mx0 <- names(mx)[which.max(mx)]
        org  <- sub("[.].*","",mx0)
        type <- sub(".*[.]","",mx0)        
        message("[probe2symbol] mapped ",format(100*max.mx,digits=2),"% of probes")
        if(max.mx < 0.5 && max.mx>0) {
            message("[probe2symbol] WARNING! low mapping ratio: r= ",max.mx)
        }        
        if(max.mx==0) {
            message("[probe2symbol] WARNING! zero mapping ratio: r= ")
            type = NULL
        }
        org
        type
    }
    if(is.null(type)) {
        cat("probe2symbol: invalid type: ",type,"\n")
        return(NULL)
    }
    if(!type %in% c("ensembl","ensemblTRANS","unigene","refseq","accnum","uniprot","symbol")) {
        cat("probe2symbol: invalid type: ",type,"\n")
        return(NULL)
    }

    cat("[probe2symbol] organism = ",org,"\n")
    cat("[probe2symbol] probe.type = ",type,"\n")
    type

    if(type=="symbol") {
        cat("probe2symbol: probe is already symbol\n")
        if(any(grep(" /// ",probes))) {
            symbol0 <- strsplit(probes, split=" /// ")
        } else if(any(grep("[;,]",probes))) {
            symbol0 <- strsplit(probes, split="[;,\\|]")
        } else {
            symbol0 <- probes
        }
        ## all.symbols <- NULL
        ## if(org=="human") all.symbols <- unlist(as.list(org.Hs.egSYMBOL))
        ## if(org=="mouse") all.symbols <- unlist(as.list(org.Mm.egSYMBOL))
        ## symbol0 <- lapply(symbol0, function(s) intersect(s,all.symbols))

    } else {
        org
        if(org=="human") {
            symbol0 <- AnnotationDbi::mapIds(org.Hs.eg.db, probes, 'SYMBOL', toupper(type))
        }
        if(org=="mouse") {
            symbol0 <- AnnotationDbi::mapIds(org.Mm.eg.db, probes, 'SYMBOL', toupper(type))
        }
    }

    ## Unrecognize probes
    nna <- which(is.na(names(symbol0)))
    length(nna)
    if(length(nna)) names(symbol0)[nna] <- probes[nna]

    ## What to do with unmapped/missing symbols????
    symbol <- sapply(symbol0,"[",1)  ## takes first symbol only!!!
    isnull <- which(sapply(symbol,is.null))
    symbol[isnull] <- NA
    if(keep.na) {
        sel.na <- which(is.na(symbol))
        symbol[sel.na] <- probes[sel.na]
    }
    symbol <- unlist(symbol)
    names(symbol) <- NULL
    Matrix::head(symbol)

    symbol    
}


##s=title

trimsame <- function(s, split=" ", ends=TRUE, summarize=FALSE) {
    if(ends) return(trimsame.ends(s, split=split, summarize=summarize))
    return(trimsame0(s, split=split, summarize=summarize))
}

trimsame.ends <- function(s, split=" ", summarize=FALSE) {
    s1 <- trimsame0(s, split=split, summarize=summarize)
    s2 <- sapply(strsplit(s1, split=split),function(x) paste(rev(x),collapse=" "))
    s3 <- trimsame0(s2, split=split, summarize=summarize, rev=TRUE)
    s4 <- sapply(strsplit(s3, split=split),function(x) paste(rev(x),collapse=" "))
    s4
}

trimsame0 <- function(s, split=" ", summarize=FALSE, rev=FALSE) {
    for(i in 1:4) s <- gsub(paste0(split,split),split,s)
    whereSpaces <- function(s) as.vector(gregexpr(split, s)[[1]])
    sp <- whereSpaces(s[1])
    sp
    if(length(sp)==0) return(s)
    
    is.same = TRUE
    i=1
    j=0
    while(is.same && i<=length(sp)) {
        is.same <- all(duplicated(substring(s,1,sp[i]))[-1])
        if(is.same) j = i
        i = i + 1
    }
    i
    j
    is.same
    s1 <- s
    if(j>0) {
        samepart <- substring(s[1],1,sp[j])
        samepart
        subst <- ""
        if(summarize) {
            subst <- substring(strsplit(trimws(samepart),split=split)[[1]],1,1)
            if(rev) subst <- rev(subst)
            subst <- paste(c(toupper(subst)," "),collapse="")
        }
        s1 <- sub(samepart,subst,s)
    }
    s1
}

##s=rep("abc",100)
dbg.BAK <- function(... ) {
    if(exists("DEBUG") && DEBUG) {
        ##msg = paste0(ifelse(is.null(module),"",paste0("<",module,"> ")),msg)
        msg = sapply( list(...),paste,collapse=" ")
        message(cat(paste0("[DBG] ",sub("\n$","",paste(msg,collapse=" ")),"\n")))
    }
}

read.csv3.BAK <- function(file, ...)
{
    ## read delimited table automatically determine separator
    line1 <- as.character(read.csv(file, comment.char='#', sep='\n',nrow=1)[1,])
    sep = names(which.max(sapply(c('\t',',',';'),function(s) length(strsplit(line1,split=s)[[1]]))))
    message("[read.csv3] sep = ",sep)
    read.csv(file, comment.char='#', sep=sep, ...)
}

##check.names=FALSE;row.names=1;stringsAsFactors=FALSE;header=TRUE
read.csv3 <- function(file, as_matrix=FALSE)
{
    ## read delimited table automatically determine separator. Avoid
    ## duplicated rownames.
    line1 <- as.character(read.csv(file, comment.char='#', sep='\n',nrow=1)[1,])
    sep = names(which.max(sapply(c('\t',',',';'),function(s) length(strsplit(line1,split=s)[[1]]))))
    ##message("[read.csv3] sep = ",sep)
    ##x <- read.csv(file, comment.char='#', sep=sep)
    sep
    ##x <- read.csv(file, comment.char='#', sep=sep, check.names=FALSE, stringsAsFactors=FALSE)
    x <- data.table::fread(file, sep=sep, check.names=FALSE, stringsAsFactors=FALSE, header=TRUE)
    x <- as.data.frame(x)
    x <- x[grep("^#",x[[1]],invert=TRUE),,drop=FALSE]  ## drop comments
    dim(x)
    xnames <- as.character(x[,1])
    sel <- which(xnames!="" & !duplicated(xnames))
    x <- x[sel,-1,drop=FALSE]
    if(as_matrix) x <- as.matrix(x)
    if(length(sel)) {
        rownames(x) <- xnames[sel]
    }
    ##x <- type.convert(x)
    x
}

read.as_matrix.SAVE <- function(file)
{
    ## read delimited table automatically determine separator. allow duplicated rownames.
    line1 <- as.character(read.csv(file, comment.char='#', sep='\n',nrow=1)[1,])
    sep = names(which.max(sapply(c('\t',',',';'),function(s) length(strsplit(line1,split=s)[[1]]))))    
    x0 <- read.csv(file, comment.char='#', sep=sep, check.names=FALSE, stringsAsFactors=FALSE)
    x <- NULL
    sel <- which(! as.character(x0[,1]) %in% c(""," ","NA","na",NA))        
    if(length(sel)) {
        x <- as.matrix(x0[sel, -1 ,drop=FALSE])  ## always as matrix
        rownames(x) <- x0[sel,1]
    }
    return(x)
}

read.as_matrix <- function(file)
{
    ## read delimited table automatically determine separator. allow
    ## duplicated rownames. This implements with faster fread.
    x0 <- data.table::fread(file=file, check.names=FALSE, header=TRUE,
                            blank.lines.skip=TRUE, stringsAsFactors=FALSE)
    x <- NULL
    sel <- which(!as.character(x0[[1]]) %in% c(""," ","NA","na",NA))
    length(sel)    
    if(length(sel)) {
        x <- as.matrix(x0[sel, -1 ,drop=FALSE])  ## always as matrix
        rownames(x) <- x0[[1]][sel]
    }
    return(x)
}

##check.names=FALSE;row.names=1;stringsAsFactors=FALSE;header=TRUE
fread.csv <- function(file, check.names=FALSE, row.names=1,
                      stringsAsFactors=FALSE, header=TRUE, asMatrix=TRUE)
{
    
    df <- data.table::fread(file=file, check.names=check.names, header=header)
    x <- data.frame(df[,2:ncol(df)], stringsAsFactors=stringsAsFactors,
                    check.names=check.names)
    is.num <- all(sapply(x,class)=="numeric")
    is.char <- all(sapply(x,class)=="character")
    is.int <- all(sapply(x,class)=="integer")
    if(asMatrix && (is.num || is.char || is.int)) x <- as.matrix(x)
    rownames(x) <- df[[row.names]]  ## allow dups if matrix
    return(x)
}

tagDuplicates <- function(s) {
    ## Tag duplicate with blanks
    ##
    jj <- which(duplicated(s))
    t <- s[jj][1]
    for(t in unique(s[jj])) {
        ii <- which(s==t)
        spaces <- paste(rep(" ",length(ii)),collapse="")
        blanks <- substring(spaces,0,0:(length(ii)-1))
        ##s[ii] <- paste(s[ii],1:length(ii),sep=".")
        s[ii] <- paste0(s[ii],blanks)
    }
    s <- gsub("[.]1$","",s)
    s
}

wrapHyperLink <- function(s, gs) {
    
    ## GEO/GSE accession
    gs = as.character(gs)
    s1 = s = as.character(s)
    jj <- grep("GSE[0-9]",gs)
    if(length(jj)) {
        acc = sub("[-_ ].*","",gsub("^.*GSE","GSE",gs[jj]))
        url = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",acc)
        s1[jj] <- paste0("<a href='",url,"' target='_blank'>",s[jj],"</a>")
    }

    ## KEGG
    jj <- grep("HSA-[0-9][0-9]|hsa[0-9][0-9]",gs)
    if(length(jj)) {
        id = sub("^.*HSA-|.*hsa","hsa",gs[jj])
        id = sub("[-_ ].*","",id)
        url = paste0("https://www.genome.jp/kegg-bin/show_pathway?map=",id,"&show_description=show")
        s1[jj] <- paste0("<a href='",url,"' target='_blank'>",s[jj],"</a>")
    }
    
    ## Reactome (afer doing KEGG???)
    jj <- grep("R-HSA-[0-9][0-9]",gs)
    if(length(jj)) {
        id = sub("^.*R-HSA-","R-HSA-",gs[jj])
        url = paste0("https://reactome.org/content/detail/",id)
        s1[jj] <- paste0("<a href='",url,"' target='_blank'>",s[jj],"</a>")
    }

    ## Wikipathways
    jj <- grep("_WP[0-9][0-9]",gs)
    if(length(jj)) {
        id = sub("^.*_WP","WP",gs[jj])
        url = paste0("https://www.wikipathways.org/index.php/Pathway:",id)
        s1[jj] <- paste0("<a href='",url,"' target='_blank'>",s[jj],"</a>")
    }

    ## MSigDB
    jj <- grep("^H:|^C[1-8]:|HALLMARK",gs)
    if(length(jj)) {
        gs1 = sub("^.*:","",gs[jj])
        url = paste0("http://software.broadinstitute.org/gsea/msigdb/cards/",gs1)
        s1[jj] <- paste0("<a href='",url,"' target='_blank'>",s[jj],"</a>")
    }
    
    ## GO reference
    jj <- grep("\\(GO_.*\\)$",gs)
    if(length(jj)) {
        id = gsub("^.*\\(GO_|\\)$","",gs[jj])
        url = paste0("http://amigo.geneontology.org/amigo/term/GO:",id)
        s1[jj] <- paste0("<a href='",url,"' target='_blank'>",s[jj],"</a>")
    }
    
    return(s1)
}

reverse.AvsB <- function(comp) {
    reverse.AvsB.1 <- function(comp) {
        prefix=postfix=""
        if(any(grepl("[:]",comp))) prefix <- sub(":.*","",comp)
        if(any(grepl("[@]",comp))) postfix <- sub(".*@","",comp)
        comp0 <- gsub(".:|@.*","",comp)
        ab <- paste(rev(strsplit(comp0, split="_vs_|_VS_")[[1]]),collapse="_vs_")
        gsub("^:|@$","",paste0(prefix,":",ab,"@",postfix))
    }
    as.character(sapply(comp,reverse.AvsB.1))
}

is.POSvsNEG <- function(pgx) {
    ## Determines automagically from contrast matrix if notation is
    ## 'A_vs_B' or 'B_vs_A' (which group is positive in the contrast
    ## matrix). Too complicated... maybe we should just require one
    ## definition...
    ##
    ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ## !!!!!!!!!!!! We should get rid of this... !!!!!!!!!!!!!
    ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    cntrmat <- pgx$model.parameters$contr.matrix
    design <- pgx$model.parameters$design
    ##ct0 <- cntrmat[,comp]        

    ## rely only on contrasts with '_vs_'
    cntrmat <- cntrmat[,grep("_vs_",colnames(cntrmat)),drop=FALSE]
    dim(cntrmat)
    grp1 <- sapply(strsplit(colnames(cntrmat),split="_vs_"),"[",1)
    grp2 <- sapply(strsplit(colnames(cntrmat),split="_vs_"),"[",2)
    grp1 <- sub(".*[:]|@.*","",grp1)
    grp2 <- sub(".*[:]|@.*","",grp2)
    is.null(design)
    is.PosvsNeg1 <- NA
    grp1
    grp2
    if(FALSE && !is.null(design)) {
        is.pn <- rep(NA,length(grp1))
        i=1
        for(i in 1:length(grp1)) {
            ##grp1x <- intersect(grp1[i],rownames(cntrmat))
            ##grp2x <- intersect(grp2[i],rownames(cntrmat))        
            ##grp1.sign <- mean(cntrmat[intersect(grp1,rownames(cntrmat)),which(grp1 %in% grp1x)])
            ##grp2.sign <- mean(cntrmat[intersect(grp2,rownames(cntrmat)),which(grp2 %in% grp2x)])
            j1 <- grep(grp1[i], rownames(cntrmat), fixed=TRUE)
            j2 <- grep(grp2[i], rownames(cntrmat), fixed=TRUE)
            grp1.sign <- mean(cntrmat[j1,i], na.rm=TRUE)
            grp2.sign <- mean(cntrmat[j2,i], na.rm=TRUE)
            grp1.sign
            grp2.sign
            if(!is.nan(grp1.sign) && !is.nan(grp2.sign)) {
                is.pn[i] <- (grp1.sign > grp2.sign )
                ##is.NegvsPos1 <- (grp2.sign > 0 && grp1.sign < 0)
            }
        }
        is.pn
        is.PosvsNeg1 <- mean(is.pn,na.rm=TRUE)>0
    } else {
        ## This uses the experiment matrix (sample-based contrast) and
        ## the annotation to determine is A_vs_B or B_vs_A was
        ## intended.
        ##
        expmat <- pgx$model.parameters$exp.matrix
        is.pn <- rep(NA,length(grp1))
        i=1
        for(i in 1:length(grp1)) {
            ##a1 <- apply(pgx$samples,1,function(a) mean(a %in% grp1[i]))
            ##a2 <- apply(pgx$samples,1,function(a) mean(a %in% grp2[i]))
            a1 <- apply(pgx$samples,1,function(a) mean(grepl(grp1[i],a)))
            a2 <- apply(pgx$samples,1,function(a) mean(grepl(grp2[i],a)))
            j1 <- which(a1 > a2)   ## samples with phenotype  more in grp1
            j2 <- which(a2 >= a1)  ## samples with phenotype  more in grp2
            s1=s2=0
            if(length(j1)) s1 <- rowMeans(expmat[j1,i,drop=FALSE] > 0, na.rm=TRUE)
            if(length(j2)) s2 <- rowMeans(expmat[j2,i,drop=FALSE] > 0, na.rm=TRUE)
            mean(s1)
            mean(s2)
            if(mean(s1) > mean(s2)) is.pn[i] <- TRUE
            if(mean(s2) > mean(s1)) is.pn[i] <- FALSE
        }
        is.pn
        is.PosvsNeg1 <- mean(is.pn,na.rm=TRUE)>0
    }
    is.PosvsNeg1
    
    ## look for keywords 
    grp1.neg2 <- mean(grepl("neg|untr|ref|wt|ctr|control",tolower(grp1)))
    grp2.neg2 <- mean(grepl("neg|untr|ref|wt|ctr|control",tolower(grp2)))
    grp1.neg2
    grp2.neg2
    is.PosvsNeg2 <- NA
    if(grp1.neg2>0 || grp2.neg2>0) {
        is.PosvsNeg2 <- (grp2.neg2 > grp1.neg2)
    }
    is.PosvsNeg2
    
    ok <- setdiff(c(is.PosvsNeg1,is.PosvsNeg2),NA)[1]  ## priority to first test??
    if(is.na(ok) || length(ok)==0) ok <- TRUE  ## DEFAULT if not known !!!!!
    ok
}


is.categorical <- function(x, max.ncat=null, min.ncat=2) {
    max.ncat <- length(x)-1
    is.factor <- any(class(x) %in% c("factor","character"))
    is.factor
    n.unique <- length(unique(setdiff(x,NA)))
    n.notna  <- length(x[!is.na(x)])
    n.unique
    n.notna    
    is.id    <- (n.unique > 0.8*n.notna)
    is.id
    is.factor2 <- (is.factor & !is.id & n.unique>=min.ncat & n.unique<= max.ncat)
    is.factor2
}

##remove.dup=TRUE;min.ncat=2;max.ncat=20
pgx.discretizePhenotypeMatrix <- function(df, min.ncat=2, max.ncat=20, remove.dup=FALSE)
{
    catpheno <- pgx.getCategoricalPhenotypes(
        df, max.ncat=max.ncat, min.ncat=min.ncat, remove.dup=remove.dup)
    numpheno <- pgx.getNumericalPhenotypes(df)
    catpheno
    numpheno
    numpheno <- setdiff(numpheno, catpheno) ## already in categories?
    df.num <- c()
    if(length(numpheno)) {
        df.num <- as.matrix(df[,numpheno,drop=FALSE])
        is.high <- t(t(df.num) > apply(df.num,2,median,na.rm=TRUE))
        df.num[is.high]  <- "high"
        df.num[!is.high] <- "low"
    }
    df1 <- df[,0]
    if(length(catpheno)) df1 <- cbind(df1, df[,catpheno,drop=FALSE])
    if(length(numpheno)) df1 <- cbind(df1, df.num)
    rownames(df1) <- rownames(df)
    df1
}

pgx.getNumericalPhenotypes <-function(df)
{
    is.bad = 0
    is.bad1 <- grepl("^sample$|[_.]id$|replic|rep|patient|donor|individ",tolower(colnames(df)))
    is.bad2 <- grepl("year|month|day|^efs|^dfs|surv|follow",tolower(colnames(df)))    
    is.bad3 <- apply(df,2,function(x) any(grepl("^sample|patient|replicate|donor|individ",x,ignore.case=TRUE)))    
    is.bad <- (is.bad1 | is.bad2 | is.bad3)
    table(is.bad)
    is.bad

    numratio <- apply(df,2,function(x)length(unique(x))) / nrow(df)
    numratio
    numpheno <- (apply(df,2,is.num) & !is.bad & numratio>0.5)
    numpheno    
    names(which(numpheno==TRUE))
}

##max.ncat=9999;min.ncat=2
pgx.getCategoricalPhenotypes <-function(df, min.ncat=2, max.ncat=20, remove.dup=FALSE) {
    ##
    ##
    ##
    ##df <- type.convert(df)
    
    is.bad = 0

    ## ... exclude sample IDs
    ##is.bad1 <- grepl("^sample$|[_.]id$|replic|rep|patient|donor|individ",tolower(colnames(df)))
    is.bad1 <- grepl("^sample$|[_.]id$|patient|donor|individ",tolower(colnames(df)))

    ## ... exclude numerical dates/age/year
    is.bad2 <- grepl("ratio|year|month|day|^age$|^efs|^dfs|surv|follow",tolower(colnames(df)))    
    ## is.factor <- sapply(sapply(df, class), function(s) any(s %in% c("factor","character")))
    is.num  <- sapply(df,class) == "numeric"
    is.bad2 <- (is.bad2 & is.num)  ## no numeric

    ## ... exclude any sample ID coded in columns...
    is.bad3 <- apply(df,2,function(x) mean(grepl("^sample|patient|donor|individ",
                                                 x[!is.na(x)],ignore.case=TRUE))>0.8)    
    ## is.bad <- (is.bad1 | is.bad2 | is.bad3)
    is.bad <- (is.bad2 | is.bad3)    
    is.bad
    table(is.bad)

    ## auto-determine which are factors
    is.factor <- apply(df, 2, is.categorical)
    is.factor
    n.unique <- apply(df,2,function(x) length(unique(setdiff(x,c(NA,"NA","")))))
    n.notna <- apply(df,2,function(x) length(x[!is.na(x)]))
    is.id <- (n.unique > 0.9*n.notna)
    is.id
    is.factor2 <- (!is.bad & is.factor & !is.id & n.unique>=min.ncat & n.unique<= max.ncat)
    is.factor2

    ## take reduced matrix
    df1 <- df[,which(is.factor2),drop=FALSE]
    dim(df1)
    nlevel <- apply(df1,2,function(x) length(unique(x)))
    nchars <- apply(df1,2,function(x) max(nchar(iconv(x, "latin1", "ASCII", sub=""))))
    df1 <- df1[,order(nlevel,-nchars),drop=FALSE]
    dim(df1)

    ##head(df1)
    if(remove.dup && ncol(df1)>1) {
        i=1
        j=2
        is.dup <- rep(FALSE,ncol(df1))
        for(i in 1:(ncol(df1)-1)) {
            is.dup[i] <- FALSE
            for(j in (i+1):ncol(df1)) {
                is.dup[i] <- is.dup[i] || all(rowSums(table(df1[,i],df1[,j])!=0)==1)
            }
        }
        is.dup
        df1 <- df1[,which(!is.dup),drop=FALSE]
    }
    colnames(df1)
}


ngs.save <- function(ngs, file, update.date=TRUE, light=TRUE, system=FALSE) {

    if(update.date||is.null(ngs$date)) ngs$date <- Sys.Date()

    if(light) {
        ## ------- make a light version
        ngs$gx.meta$outputs <- NULL
        ngs$gset.meta$outputs <- NULL
        ngs$model.parameters$efit <- NULL
        ngs$gmt.all <- NULL
        ngs$families <- NULL
        ngs$collections <- NULL
        ## ngs$counts <- NULL
        ngs$gset.meta$matrices <- NULL
    }
    if(system==FALSE) {
        ## remove system (is big...)
        ngs$omicsnet <- NULL
        ngs$omicsnet.reduced <- NULL
    }
    sort(sapply(ngs, object.size)) / 1e9
    sum(sapply(ngs, object.size)) / 1e9        
    
    cat(">>> saving PGX file to",file,"\n")
    file <- iconv(file, from = '', to = 'ASCII//TRANSLIT')
    save(ngs, file=file)
}

##comp=1;level="geneset";probe=rownames(ngs$gsetX)[1]
getLevels <- function(Y) {
    yy = Y[,grep("title|name|sample|patient",colnames(Y),invert=TRUE),drop=FALSE]   ## NEED RETHINK!!!!
    is.grpvar  = apply(yy,2,function(y) max(table(y))>1)
    is.numeric = apply(yy,2,function(y) (length(table(y))/length(y)) > 0.5)
    is.grpvar = is.grpvar & !is.numeric
    yy = yy[,is.grpvar,drop=FALSE]
    ##yy = yy[,which(apply(yy,2,function(x) length(unique(x)))<20),drop=FALSE]
    ##levels = setdiff(unique(as.vector(apply(yy,2,as.character))),c("",NA))
    levels = lapply(1:ncol(yy), function(i) unique(paste0(colnames(yy)[i],"=",yy[,i])))
    levels = sort(unlist(levels))
    return(levels)
}

selectSamplesFromSelectedLevels <- function(Y, levels)
{
    if( is.null(levels) || levels=="") {
        return(rownames(Y))
    }
    pheno <- sapply(strsplit(levels,split="="),"[[",1)
    ptype <- sapply(strsplit(levels,split="="),"[[",2)
    sel = rep(1, nrow(Y))
    i=1
    for(ph in unique(pheno)) {
        k <- which(pheno==ph)
        sel <- sel * (Y[,ph] %in% ptype[k])
    }
    rownames(Y)[which(sel==1)]
}

##fields=c("symbol","name","alias","map_location","summary")
getMyGeneInfo <- function(eg, fields=c("symbol","name","alias","map_location","summary"))
{
    
    ##res = biomaRt::getGene(eg, fields="all")[[1]]
    ##fields = c("symbol","name","alias","map_location","summary")
    info = lapply(fields, function(f) biomaRt::getGene(eg, fields=f)[[1]] )
    names(info) <- fields
    info = lapply(info, function(x) ifelse(length(x)==3,x[[3]],"(not available)") )
    info = sapply(info, paste, collapse=",")
    if(0) {
        rifs = biomaRt::getGene(eg, fields="generif")[[1]]
        collapse.rif <- function(r) paste0(r$text," (PMID=",r$pubmed,")")
        rifs = rifs[which(sapply(rifs,length)>2)]
        xrif = sapply(rifs[[3]], function(rr) collapse.rif(rr))
    }
    return(info)
}

## much faster and off-line
getHSGeneInfo <- function(eg, as.link=TRUE) {
         
    env.list <- c("symbol" = org.Hs.eg.db::org.Hs.egSYMBOL,
                  "name" = org.Hs.eg.db::org.Hs.egGENENAME,
                  "map_location" = org.Hs.eg.db::org.Hs.egMAP,
                  "OMIM" = org.Hs.eg.db::org.Hs.egOMIM,
                  "KEGG" = org.Hs.eg.db::org.Hs.egPATH,
                  ##"PMID" = org.Hs.eg.db::org.Hs.egPMID,
                  "GO" = org.Hs.eg.db::org.Hs.egGO)

    info <- lapply(env.list, function(env) mget(eg, envir=env, ifnotfound=NA)[[1]])
    names(info) <- names(env.list)
    gene.symbol <- toupper(mget(as.character(eg), envir=org.Hs.eg.db::org.Hs.egSYMBOL))[1]
    info[["symbol"]] <- gene.symbol
    
    ## create link to GeneCards
    if(as.link) {
        genecards.link = "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=GENE' target='_blank'>GENE</a>"
        info[["symbol"]] <- gsub("GENE", info[["symbol"]], genecards.link)
    }

    ## create link to OMIM
    if(as.link) {
        omim.link = "<a href='https://www.omim.org/entry/OMIM' target='_blank'>OMIM</a>"
        info[["OMIM"]] <- sapply(info[["OMIM"]], function(x) gsub("OMIM",x,omim.link))
    }

    ## create link to KEGG
    kegg.link = "<a href='https://www.genome.jp/kegg-bin/show_pathway?map=hsaKEGGID&show_description=show' target='_blank'>KEGGNAME (KEGGID)</a>"
    for(i in 1:length(info[["KEGG"]])) {
        kegg.id = info[["KEGG"]][[i]]
        kegg.id = setdiff(kegg.id,NA)
        if(length(kegg.id)>0) {
            kegg.name = mget(kegg.id, envir=KEGG.db::KEGGPATHID2NAME, ifnotfound=NA)[[1]]
            if(!is.na(kegg.name) && as.link) {
                info[["KEGG"]][[i]] <- gsub("KEGGNAME",kegg.name,gsub("KEGGID",kegg.id,kegg.link))
            } else {
                info[["KEGG"]][[i]] <- kegg.name
            }
        }
    }

    ## create link to GO
    if(!is.na(info[["GO"]][1])) {
        go.evidence = c("EXP","IDA","IPI","IMP","IGI","IEP")
        amigo.link = "<a href='http://amigo.geneontology.org/amigo/term/GOID' target='_blank'>GOTERM (GOID)</a>"
        sel <- which( sapply(info[["GO"]],"[[",2) %in% go.evidence  &
                      sapply(info[["GO"]],"[[",3) %in% c("BP")  )
        sel
        info[["GO"]] <- info[["GO"]][sel]

        ## sometimes GO.db is broken...
        suppressWarnings( try.out <- try(AnnotationDbi::Term(mget("GO:0000001", envir=GO.db::GOTERM,
                                                                  ifnotfound=NA)[[1]])))
        go.ok <- (class(try.out) !="try-error")
        if(go.ok && length(sel)>0) {
            i=1
            for(i in 1:length(info[["GO"]])) {
                go_id = info[["GO"]][[i]][[1]]
                go_term = AnnotationDbi::Term(mget(go_id, envir=GO.db::GOTERM, ifnotfound=NA)[[1]])
                if(as.link) {
                    info[["GO"]][[i]] = gsub("GOTERM",go_term,gsub("GOID",go_id,amigo.link))
                } else {
                    info[["GO"]][[i]] = go_term
                }
            }
        } else {
            info[["GO"]] <- NULL
        }
    }

    return(info)
}

##levels="gene";contrast="Bmem_activation";layout=NULL;gene="IRF4";layout="layout_with_fr";hilight=NULL
pgx.getGeneFamilies <- function(genes, FILES="../files", min.size=10, max.size=500)
{

    ##dir="/home/share/datasets/gmt/";nrows=-1
    read.gmt <- function(file, dir=NULL, add.source=FALSE, nrows=-1) {
        f0 <- file
        if(strtrim(file,1)=="/") dir=NULL
        if(!is.null(dir)) f0 <- paste(sub("/$","",dir),"/",file,sep="")
        ##cat("reading GMT from file",file,"\n")
        gmt <- read.csv(f0,sep="!",header=FALSE,comment.char="#",nrows=nrows)[,1]
        gmt <- as.character(gmt)
        gmt <- gsub("[\t]+","\t",gmt)
        gmt <- sapply(gmt,strsplit,split="\t")
        names(gmt) <- NULL
        gmt.name <- sapply(gmt,"[",1)
        gmt.source <- sapply(gmt,"[",2)
        gmt.genes <- sapply(gmt,function(x) paste(x[3:length(x)],collapse=" "))
        ##gmt.genes <- gsub("[\t]+"," ",gmt.genes)
        gset <- sapply(gmt.genes,strsplit,split=" ")
        gset <- lapply(gset, function(x) setdiff(x,c("","NA",NA)))
        names(gset) <- gmt.name
        if(add.source) {
            names(gset) <- paste0(names(gset)," (",gmt.source,")")
        }
        gset <- gset[which(lapply(gset,length)>0)]
        gset
    }

    ##-----------------------------------------------------------------------------
    ## Gene families
    ##-----------------------------------------------------------------------------
    ##xgene = genes[,"gene_name"]
    ##xtitle = genes[,"gene_title"]
    families <- list()
    families[["<all>"]] = genes  ## X is sorted

    gmt.kea  <- read.gmt(file.path(FILES,"gmt/kinase_substrates_kea.gmt"))
    gmt.chea <- read.gmt(file.path(FILES,"gmt/tf_targets_chea.gmt"))
    families[["Kinases (KEA)"]] = names(gmt.kea)
    families[["Transcription factors (ChEA)"]] = names(gmt.chea)

    ## Read standard HGNC gene families (www.genefamilies.org)
    ##gmt.hgnc <- read.gmt(file.path(FILES,"gmt/hgnc_genefamilies.gmt"))
    gmt.hgnc <- read.gmt(file.path(FILES,"gmt/hgnc_genefamilies_EDITED.gmt"))

    gmt.hgnc.size <- sapply(gmt.hgnc,length)
    gmt.hgnc <- gmt.hgnc[ which( gmt.hgnc.size >=50 & gmt.hgnc.size <= 1000)]
    length(gmt.hgnc)
    names(gmt.hgnc) <- paste0(names(gmt.hgnc)," (HGNC)")

    families <- c(families, gmt.hgnc)

    ##xgene = as.character(genes$gene_name)
    ##gtype = as.character(genes$gene_type)
    ##gtype = gsub("TR_.*gene","TR_gene",gtype)
    ##gtype = gsub(".*pseudogene","pseudogene",gtype)
    ##sort(table(gtype))
    ##families[["<protein_coding>"]] = genes[which(gtype=="protein_coding")]
    ##families[["Antisense genes"]] = genes[which(gtype=="antisense")]
    ##families[["Pseudogenes"]] = genes[which(gtype=="pseudogene")]
    ##families[["TR genes"]] = genes[which(gtype=="TR_gene")]
    ##families[["lincRNA"]] = genes[which(gtype=="lincRNA")]
    ##families[["snoRNA"]] = genes[which(gtype=="snoRNA")]
    ##families[["Arginine related"]] = genes[grep("arginine",tolower(xtitle))]
    ##families[["CD family"]] = genes[grep("^CD[1-9]",genes)]
    families[["Interleukins (IL)"]] = genes[grep("^IL[1-9]",genes)]
    families[["Chemokines"]] = genes[grep("CCL|CCR|CXCR|CXCL|XCL|CX3",genes)]
    families[["Ribosomal proteins"]] = genes[grep("^RPS|^RPL",genes)]
    families[["Ribosomal (mitochondrial)"]] = genes[grep("^MRPL|^MRPS",genes)]
    ## families[["TR family"]] = genes[grep("^TR[ABDG][VJC]",genes)]
    families[["G-protein family"]] = genes[grep("^GN[ABG]|^GPBAR|^GPER|^FZD",genes)]
    families[["Heatshock proteins"]] = genes[grep("^HSP",genes)]
    families[["Integrin family"]] = genes[grep("^ITG",genes)]
    families[["MAPK family"]] = genes[grep("^MAP[1-9]|^MAPK",genes)]
    ##families[["Olfactory receptors"]] = genes[grep("^OR[1-9]",genes)]
    families[["Myosins"]] = genes[grep("^MYO[1-9]",genes)]
    families[["Protein phosphatase (PPP)"]] = genes[grep("^PPP",genes)]
    families[["PTP family"]] = genes[grep("^PTP",genes)]
    families[["Small nucleolar RNA"]] = genes[grep("^SNOR",genes)]
    families[["Toll-like receptors"]] = genes[grep("^TLR",genes)]
    ##families[["Zinc-fingers C2H2-type"]] = genes[grep("^ZNF",genes)]
    families[["EIF factors"]] = genes[grep("^EIF",genes)]
    families[["GPR proteins"]] = genes[grep("^GPR",genes)]
    families[["CDCC proteins"]] = genes[grep("^CDCC",genes)]
    families[["KIAA proteins"]] = genes[grep("^KIAA",genes)]
    families[["TRIM proteins"]] = genes[grep("^TRIM",genes)]
    families[["TNF proteins"]] = genes[grep("^TNF",genes)]
    families[["SLC proteins"]] = genes[grep("^SLC",genes)]
    families[["ZBTB proteins"]] = genes[grep("^ZBTB",genes)]
    families[["TMEM family"]] = genes[grep("^TMEM",genes)]
    families[["STAT family"]] = genes[grep("^STAT",genes)]
    families[["DNA/RNA polymerases"]] = genes[grep("^POL",genes)]
    families[["Proteasome"]] = genes[grep("^PSM",genes)]
    families[["IFN/IFIT family"]] = genes[grep("^IFN|^IFIT",genes)]
    families[["Nuclear receptors"]] = genes[grep("^NR[0-9]|^RXR|^ESR|^PGR$|^AR$|^HNF4|^ROR|^PPAR|^THR|^VDR", genes)]
    families[["Cytochrome family"]] = genes[grep("^CYP|^CYB|^CYC|^COX|^COA",genes)]
    families[["Micro RNA"]] = genes[grep("^MIR",genes)]

    ## add pathways?
    ##kk = grep("^BIOCARTA_",names(ngs$gmt.all))
    ##kk = Matrix::head(kk,5) ## just try
    ##pathways = ngs$gmt.all[kk]
    ##families = c(families, pathways)

    ## convert to mouse???
    is.mouse = (mean(grepl("[a-z]",sub(".*:","",genes))) > 0.8)
    is.mouse
    if(is.mouse) {
        require(org.Mm.eg.db)
        mouse.genes = as.character(unlist(as.list(org.Mm.egSYMBOL)))
        names(mouse.genes) = toupper(mouse.genes)
        families <- parallel::mclapply(families, function(s)
            setdiff(as.character(mouse.genes[toupper(s)]),NA),
            mc.cores=16)
    }

    ## sort all
    families = lapply(families, function(f) intersect(f,genes))

    ## ----------- filter on size
    nsize <- sapply(families, length)
    sel <- which(nsize >= min.size & nsize < max.size)
    sel <- unique(c(which(names(families)=="<all>"),sel))
    families <- families[sel]

    return(families)
}

pgx.getGeneSetCollections <- function(gsets, min.size=10, max.size=500)
{
    ##-----------------------------------------------------------------------------
    ## Gene set collections
    ##-----------------------------------------------------------------------------
    ##gsets = rownames(ngs$gsetX)

    kegg0 <- sort(gsets[grep("KEGG|.*hsa[0-9]{5}$",gsets)])
    kegg0 <- kegg0[!duplicated(getKeggID(kegg0))]

    collections = list(
        "Hallmark collection" = gsets[grep("HALLMARK",gsets)],
        "KEGG pathways" = kegg0,
        "KEGG metabolic pathways" = kegg0[grep("^00|^01",getKeggID(kegg0))],
        ##"REACTOME pathways" = gsets[grep("^REACTOME",gsets)],
        ##"Gene Ontology" = gsets[grep("^GO[_:]",gsets)],
        ##"GO biological processes" = gsets[grep("^GOBP",gsets)],
        ##"GO cellular compartment" = gsets[grep("^GOCC",gsets)],
        ##"GO molecular function" = gsets[grep("^GOMF",gsets)],
        ##"MSigDB C2" = gsets[grep("^C2",gsets)],
        ##"MSigDB C3" = gsets[grep("^C3",gsets)],
        ##"MSigDB C6" = gsets[grep("^C6",gsets)],
        ##"MSigDB C7" = gsets[grep("^C7",gsets)],
        ##"Kinase substrates (KEA)" = gsets[grep("^KEA",gsets)],
        ##"TF targets (ChEA)" = gsets[grep("^CHEA",gsets)],
        ## "Drug signatures" = gsets[grep("DRUG",gsets,ignore.case=TRUE)],
        ##"Drug signatures (up)" = gsets[grep("^DRUG_.*_UP$",gsets)],
        ##"Drug signatures (down)" = gsets[grep("^DRUG_.*_DN$",gsets)],
        ##"PID pathways" = gsets[grep("^PID_",gsets)],
        ##"PANTHER pathways" = gsets[grep("PANTHER",gsets)],
        ##"BIOCARTA pathways" = gsets[grep("^BIOCARTA_",gsets)],
        "Pathway related" = gsets[grep("pathway",gsets,ignore.case=TRUE)],
        ##"Kinase related" = gsets[grep("kinase",gsets,ignore.case=TRUE)],
        "Metabolism related" = gsets[grep("metaboli",gsets,ignore.case=TRUE)],
        "Signalling related" = gsets[grep("signal",gsets,ignore.case=TRUE)],
        "T-cell related" = gsets[grep("tcell|t-cell|t[ ]cell",gsets,ignore.case=TRUE)],
        "B-cell related" = gsets[grep("bcell]b-cell|b[ ]cell",gsets,ignore.case=TRUE)],
        "Response related" = gsets[grep("response",gsets,ignore.case=TRUE)],
        "Cancer related" = gsets[grep("cancer",gsets,ignore.case=TRUE)],
        "Immune related" = gsets[grep("immune",gsets,ignore.case=TRUE)],
        "Cell differentiation" = gsets[grep("differentiation",gsets,ignore.case=TRUE)],
        "Checkpoint related" = gsets[grep("checkpoint",gsets,ignore.case=TRUE)],
        "IL gene sets" = gsets[grep("IL[1-9]{1,2}",gsets,ignore.case=TRUE)]
    )

    collections[["<all>"]] = gsets  ## X is sorted
    collections = collections[which(sapply(collections,length)>=10)]
    collections = collections[order(names(collections))]

    ## ----------- add main collections from gene set prefixes
    gsets.db = sub(":.*","", gsets)
    gsets.groups = tapply(gsets, gsets.db, list)
    collections <- c(collections, gsets.groups)

    ## ----------- filter on size
    nsize <- sapply(collections, length)
    sel <- which( nsize >= min.size & nsize < max.size )
    collections <- collections[sel]
    return(collections)
}


##-----------------------------------------------------------------------------
## Generic module functions
##-----------------------------------------------------------------------------

filterFamily <- function(genes, family, ngs) {
    ##ngs <- shiny::isolate(inputData())
    gg = ngs$families[[10]]
    gg = ngs$families[[family]]
    ## check probe name, short probe name or gene name for match
    p0 = (toupper(sub(".*:","",rownames(genes))) %in% toupper(gg))
    p1 = (toupper(rownames(genes))  %in% toupper(gg))
    p2 = (toupper(as.character(genes$gene_name))  %in% toupper(gg))
    jj = which(p0 | p1 | p2 )
    rownames(genes)[jj]
}

filterProbes <- function(genes, gg) {
    ##genes = ngs$families[[10]]
    ##genes = ngs$families[[family]]
    ## check probe name, short probe name or gene name for match
    p0 = (toupper(sub(".*:","",rownames(genes))) %in% toupper(gg))
    p1 = (toupper(rownames(genes))  %in% toupper(gg))
    p2 = (toupper(as.character(genes$gene_name))  %in% toupper(gg))
    jj = which(p0 | p1 | p2)
    if(length(jj)==0) return(NULL)
    return( rownames(genes)[jj] )
}

computeFeatureScore <- function(X, Y, features)
{
    ##features=X=NULL
    ##Y = ngs$Y[,grep("group|sample|patient",colnames(ngs$Y),invert=TRUE)]
    sdx = apply(X,1,sd)
    names(sdx) = rownames(X)
    S = matrix(NA, nrow=length(features), ncol=ncol(Y))
    rownames(S) = names(features)
    colnames(S) = colnames(Y)
    for(k in 1:ncol(Y)) {
        grp = Y[colnames(X),k]
        grp = as.character(grp)
        score = rep(NA, length(features))
        names(score) = names(features)
        i=1
        for(i in 1:length(features)) {
            pp = features[[i]]
            ## if(input$cl_level=="gene") pp = filterProbes(ngs$GENES, features[[i]])
            pp = Matrix::head(pp[order(-sdx[pp])],100)
            mx = t(apply(X[pp,], 1, function(x) tapply(x,grp,mean)))
            ##D = as.matrix(dist(t(mx)))
            D = 1 - stats::cor(mx, use="pairwise")
            diag(D) = NA
            score[i] = mean(D,na.rm=TRUE)
        }
        S[,k] = score
    }
    if(is.null(S)) return(NULL)
    return(S)
}


##-----------------------------------------------------------------------------
## KEGG
##-----------------------------------------------------------------------------

getKeggID <- function(gsets)
{
    ## Guess KEGG id from gene set name
    ##

    kegg.names <- unlist(as.list(KEGG.db::KEGGPATHID2NAME))
    kegg.ids <- names(kegg.names)
    kegg.namesUPPERCASE <- toupper(gsub("[- ]","_",gsub("[)(/,.']","",kegg.names)))
    kegg.namesUPPERCASE <- gsub("__|___","_",kegg.namesUPPERCASE)

    k1 <- grep("hsa[0-9]{5}$",gsets,value=TRUE,ignore.case=TRUE)
    names(k1) <- sub(".*_hsa0","0",k1)
    k2 <- grep("kegg",gsets,value=TRUE,ignore.case=TRUE)
    k2x <- sub(".*KEGG_","",k2)
    names(k2) <- kegg.ids[match(k2x,kegg.namesUPPERCASE)]

    kegg.gsets <- c(k1, k2)
    kegg.ids <- names(kegg.gsets)[match(gsets, kegg.gsets)]
    kegg.ids

    return(kegg.ids)
}

##-----------------------------------------------------------------------------
## Generic helper functions
##-----------------------------------------------------------------------------

is.Date <- function(x) {
    ## From https://stackoverflow.com/questions/18178451/is-there-a-way-to-check-if-a-column-is-a-date-in-r
    if (!all(is.na(as.Date(
             as.character(x),
             format = c("%d/%m/%Y", "%d-%m-%Y", "%Y/%m/%d", "%Y-%m-%d")
         )))) {
        return(TRUE)
    } else{
        return(FALSE)
    }
}

##mat=ngs$X;group=ngs$samples$group;FUN=mean
averageByGroup <- function(mat, group, FUN=mean) {
    ##sum(duplicated(ngs$genes$gene_name))
    ##out <- t(apply(mat, 1, function(x) tapply(x, group, FUN)))
    out <- do.call(cbind, tapply(1:ncol(mat), group,
                                 function(i) rowMeans(mat[,i,drop=FALSE])))
    return(out)
}

makeAcronym <- function(x) {
    xp <- strsplit(x,split="[_ -]")
    sapply(xp, function(s) {
        if(length(s)==1) return(substring(s,1,2))
        toupper(paste(substring(s,1,1),collapse=""))
    })
}


relevelFactorFirst <- function(f) {
    factor(f, levels=f[!duplicated(f)])
}

extremeCorrelation <- function(query_sig, ref_set, n=200) {
    gg <- intersect(rownames(ref_set), names(query_sig))
    if(n>0) {
        gg <- gg[unique(c(Matrix::head(order(query_sig),n),head(order(-query_sig),n)))]
    }
    rho <- stats::cor( ref_set[gg,], query_sig[gg], use="pairwise")
    rho <- rho[order(-rowMeans(rho**2,na.rm=TRUE)),,drop=FALSE]
    if(NCOL(rho)==1) rho <- rho[,1]
    return(rho)
}

##s=symbol;org="hs"
alias2hugo <- function(s, org=NULL, na.orig=TRUE) {
    
    require(org.Hs.eg.db)
    require(org.Mm.eg.db)

    hs.symbol <- unlist(as.list(org.Hs.egSYMBOL))
    mm.symbol <- unlist(as.list(org.Mm.egSYMBOL))
    if(is.null(org)) {
        is.human <- mean(s %in% hs.symbol,na.rm=TRUE) > mean(s %in% mm.symbol,na.rm=TRUE)
        org <- ifelse(is.human,"hs","mm")
    }
    org    
    ##eg <- sapply(lapply(s, get, env=org.Hs.eg.db::org.Hs.egALIAS2EG),"[",1)
    nna = which(!is.na(s) & s!="" & s!=" ")
    s1 <- trimws(s[nna])
    hugo <- NULL
    if(org == "hs") {
        eg <- sapply(mget(s1, envir=org.Hs.egALIAS2EG, ifnotfound=NA),"[",1)
        eg[is.na(eg)] <- "unknown"
        hugo <- sapply(mget(eg, envir=org.Hs.egSYMBOL, ifnotfound=NA),"[",1)
    } else if(org == "mm") {
        eg <- sapply(mget(s1, envir=org.Mm.egALIAS2EG, ifnotfound=NA),"[",1)
        eg[is.na(eg)] <- "unknown"
        hugo <- sapply(mget(eg, envir=org.Mm.egSYMBOL, ifnotfound=NA),"[",1)        
    } else {
        stop("[alias2hugo] invalid organism")
    }
    jj <- which(is.na(hugo))
    if(na.orig && length(jj)) hugo[jj] <- s1[jj]
    hugo0 <- rep(NA,length(s))
    hugo0[nna] <- hugo
    return(hugo0)
}

breakstringBROKEN <- function(s, n, force=FALSE) {
    if(is.null(s) || length(s)==0) return(NULL)
    if(length(s)==1) return(breakstring1(s,n=n,force=force))
    sapply(s, breakstring1,n=n,force=force)
}

##s="breakstringBROKENbreakstringBROKENbreakstringBROKENbreakstringBROKEN";n=10
breakstring <- function(s, n, nmax=999, force=FALSE, brk='\n') {
    if(is.na(s)) return(NA)
    s <- substring(as.character(s),1,nmax)
    if(nchar(s)<n) return(s)
    b = substring(s,1,n)
    n1 <- n+1
    for(i in 1:10) {
        if(n1>nchar(s)) break
        b1 <- substring(s,n1,n1+n)
        b <- paste0(b,brk,b1)
        n1 <- n1+n+1
    }
    return(b)
}

##s="breakstringBROKENbreakstringBROKENbreakstringBROKENbreakstringBROKEN";n=10
##n=20;brk="\n"
breakstring2 <- function(s, n, brk='\n', nmax=999) {
    if(is.na(s)) return(NA)
    s <- substring(as.character(s),1,nmax)
    if(is.na(s)) return(NA)
    if(nchar(s)<n) return(s)
    a = ""
    words <- paste0(strsplit(s, split=" ")[[1]]," ")
    words
    it =0
    len1 = sum(sapply(words,nchar))
    len1
    while(len1>n && it<100 && length(words)>1) {
        len = cumsum(sapply(words,nchar))
        len
        k = which(len>n)[1]
        k
        if(k==1) k=2
        a <- paste0(a,paste0(words[1:(k-1)],collapse=""),brk)
        a
        words <- words[k:length(words)]
        words
        it = it+1
        len1 = sum(sapply(words,nchar))
        len1
    }
    a <- paste0(a,paste0(words,collapse=""),brk)
    a <- sub(paste0(brk,"$"),"",gsub(paste0(" ",brk),brk,a))
    return(a)
}

shortstring <- function(s, n, dots=1 ) {
    sapply(s, shortstring0, n=n, dots=dots)
}

shortstring0 <- function(s, n, dots=1 ) {
    ##s0 = as.character(s)
    s0 <- iconv(as.character(s),to="UTF-8")
    s0 <- gsub("[&].*[;]","",s0) ## HTML special garbage...
    jj <- which(nchar(s0)>n)
    if(length(jj)==0) return(s0)
    s <- s0[jj]
    if(dots<1) {
        n1 <- ceiling(dots*n)
        s1 <- substring(s,1, n1 )
        s2 <- substring(s, nchar(s) - n + nchar(s1),nchar(s))
        aa <- paste0(s1,"...",s2)
    } else {
        aa <- paste0(substring(s,1,n),"...")
    }
    s0[jj] <- aa
    s0
}

psort <- function(x,p.col=NULL) {
    j = grep("p.value|^p$|p-val|pval",tolower(colnames(x)))[1]
    x[order(x[,j]),]
}

color_from_middle <- function (data, color1,color2) {
    ## from https://stackoverflow.com/questions/33521828/
    max_val=max(abs(data),na.rm=TRUE)
    DT::JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",max_val,color1,max_val,color1,color2,color2,max_val,max_val))
}

tidy.dataframe <- function(Y) {
    
    ##as_tibble(Y)
    Y <- Y[,which(colMeans(is.na(Y))<1),drop=FALSE]
    Y <- apply(Y,2,function(x) sub("^NA$",NA,x)) ## all characters
    Y <- Y[,which(colMeans(is.na(Y))<1),drop=FALSE]
    Y <- apply(Y,2,function(x) gsub("^[ ]*|[ ]*$","",x))
    suppressWarnings( num.Y <- apply(Y,2,function(x) as.numeric(as.character(x))) )
    is.numeric <- ( 0.8*colMeans(is.na(num.Y)) <= colMeans(is.na(Y)) )
    nlevel <- apply(Y,2,function(x) length(unique(x)))
    is.factor  <- (!is.numeric | (is.numeric & nlevel<=3) )
    is.factor <- ( is.factor | grepl("batch|replicat|type|clust|group",colnames(Y)) )
    new.Y <- data.frame(Y, check.names=FALSE)
    new.Y[,which(is.numeric)] <- num.Y[,which(is.numeric),drop=FALSE]
    new.Y[,which(is.factor)]  <- apply(Y[,which(is.factor),drop=FALSE],2,
                                       function(a) factor(as.character(a)))
    new.Y <- data.frame(new.Y, check.names=FALSE)
    return(new.Y)
}

param.class <- function(A) sapply(tidy.dataframe(A),class)

is.num <- function(y) {
    suppressWarnings(numy <- as.numeric(as.character(y)))
    t1 <- !all(is.na(numy)) && is.numeric(numy)
    t2 <- length(unique(y)) > 0.33*length(y)
    (t1 && t2)
}

isanumber <- function(x) {
    x <- sub("NA",NA,x)
    x[which(x=="")] <- NA
    suppressWarnings(nx <- as.numeric(x[!is.na(x)]))
    (length(nx)>0 && mean(!is.na(nx)) > 0.5)
}

expandAnnotationMatrix <- function(A) {
    expandPhenoMatrix(A)
}

expandAnnotationMatrixSAVE <- function(A) {
    ## get expanded annotation matrix
    nlevel <- apply(A,2,function(x) length(unique(x)))
    y.isnum <- apply(A,2,is.num)
    ##kk <- (y.isnum | (!y.isnum & nlevel>1 & nlevel<=5))
    ##A <- A[,which(kk),drop=FALSE]
    Matrix::head(A)
    i=1
    m1 <- list()
    for(i in 1:ncol(A)) {
        if(is.num(A[,i])) {
            m0 <- matrix(rank(A[,i],na.last="keep"), ncol=1)
            colnames(m0) <- colnames(A)[i]
        } else {
            x <- as.character(A[,i])
            x[is.na(x)] <- "_"
            m0 <- model.matrix( ~ 0 + x)
            colnames(m0) <- sub("^x","",colnames(m0))
        }
        if(NCOL(m0)>1) {
            colnames(m0) <- paste0(colnames(A)[i],"=",colnames(m0))
        }
        m1[[i]] <- m0
    }
    names(m1) <- colnames(A)
    ##m1 <- lapply(m1, function(m) {colnames(m)=sub("^x","",colnames(m));m})
    M <- do.call( cbind, m1)
    rownames(M) <- rownames(A)
    return(M)
}

expandPhenoMatrix <- function(pheno, collapse=TRUE, drop.ref=TRUE) {
    ## get expanded annotation matrix
    ##a1 <- pheno[,grep("group|sample",colnames(pheno),invert=TRUE),drop=FALSE]
    ##m2 <- expandAnnotationMatrix(a1)
    a1 <- tidy.dataframe(pheno)
    nlevel <- apply(a1,2,function(x) length(setdiff(unique(x),NA)))
    nterms <- colSums(!is.na(a1))
    ##y.class <- sapply(a1,class)
    y.class <- sapply(type.convert(pheno),class)
    ##y.isnum <- apply(a1,2,is.num)
    y.isnum <- (y.class %in% c("numeric","integer"))
    nlevel
    nratio <- nlevel/nterms
    nratio
    kk <- which(y.isnum | (!y.isnum & nlevel>1 & nratio < 0.66 ))
    if(length(kk)==0) return(NULL)
    a1 <- a1[,kk,drop=FALSE]
    a1.isnum <- y.isnum[kk]
    Matrix::head(a1)
    i=1
    m1 <- list()
    for(i in 1:ncol(a1)) {
        if( a1.isnum[i] ) {
            ##m0 <- matrix(rank(a1[,i]), ncol=1)
            suppressWarnings( x <- as.numeric(a1[,i]) )
            m0 <- matrix( (x > median(x,na.rm=TRUE)), ncol=1)
            colnames(m0) <- "high"
        } else if(drop.ref && nlevel[i]==2) {
            x <- as.character(a1[,i])
            x1 <- Matrix::tail(sort(x),1)
            m0 <- matrix(x==x1, ncol=1)
            colnames(m0) <- x1
        } else {
            x <- as.character(a1[,i])
            x[is.na(x) | x=="NA" | x==" "] <- "_"
            m0 <- model.matrix( ~ 0 + x)
            colnames(m0) <- sub("^x","",colnames(m0))
        }
        rownames(m0) <- rownames(a1)
        ## remove "_"
        if("_" %in% colnames(m0)) {
            m0 <- m0[,-which(colnames(m0)=="_")]
        }

        m1[[i]] <- m0
    }
    names(m1) <- colnames(a1)
    if(collapse) {
        for(i in 1:length(m1)) {
            colnames(m1[[i]]) <- paste0(names(m1)[i],"=",colnames(m1[[i]]))
        }
        m1 <- do.call(cbind, m1)
    }
    return(m1)
}

correctMarchSeptemberGenes <- function(gg) {
    sep.from <- c( paste0("0",1:9,"-Sep"), paste0(1:19,"-Sep"))
    sep.to <- c( paste0("SEPT",1:9), paste0("SEPT",1:19))
    mar.from <- c( paste0("0",1:9,"-Mar"), paste0(1:19,"-Mar"))
    mar.to <- c( paste0("MARCH",1:9), paste0("MARCH",1:19))

    
    from <- c(sep.from, mar.from)
    to <- c(sep.to, mar.to)
    gg1 <- sub("[.-]Sep$|[.-]SEP$","-Sep",gg)
    gg1 <- sub("[.-]Mar$|[.-]MAR$","-Mar",gg1)
    jj <- which( from %in% gg1)
    gg2 <- gg1
    if(length(jj)>0) {
        cat("Found ",length(jj),"Sept/Mar genes!\n")
        length(jj)
        from[jj]
        gg2 <- plyr::mapvalues(gg1, from[jj], to[jj])
    }
    return(gg2)
}

cor.pvalue <- function(x,n) pnorm(-abs(x/((1-x**2)/(n-2))**0.5))

##=====================================================================================
##=========================== END OF FILE =============================================
##=====================================================================================
