library(plotly)
library(manhattanly)



####
#VOLCANO
###
#A <-   matrix( c(2, 4, 14, 1, 5, 16, 0, 12, 14, 3, 4, 5,1:12),  nrow=3, ncol=6)
#rownames(A) <- c("BAZ1B", "GAPDH", "ENO1")
##a=P[, group1];b=P[, group2];c=genes
volcano <- function(a, b, c){
    ## a is the first group
    ## b is the second group
    ## c is rownames of A
    B <- cbind(a,b)
    pvals = apply(B,1,function(x) {t.test(x[1:ncol(a)],x[(ncol(a)+1):(ncol(a)+ncol(b))])$p.value})
    if(0) {
        pvals <- rep(NA,nrow(B))
        for(i in 1:nrow(B)) {
            x = B[i,]
            pv <- t.test(x[1:ncol(a)],x[(ncol(a)+1):(ncol(a)+ncol(b))])$p.value
            pvals[i] <- pv
        }
        ##pvals = apply(B,1,function(x) {t.test(x[1:ncol(a)],x[(ncol(a)+1):(ncol(a)+ncol(b))])$p.value})






    }

    diffs = apply(a,1,mean, na.rm=T)-apply(b,1,mean,na.rm=T)
    data = data.frame(cbind(diffs,pvals))
    data$Gene <- c
    colnames(data) <- c("EFFECTSIZE", "P", "Gene")
    return(data)
}

#data <- volcano(A[,1:3],A[,4:6], rownames(A))
#volcanoly(data, gene = "Gene")

###
# COPY NUMBER ESTIMATION
###

CopyNumber <- function(data, mol,y){
	# data should be a numeric matrix
	# mol is the molecular weight
	# y is the mass of the cell in PICOGRAM (pg)
	TotalIntenisty <- apply(data, 2, sum)
	# mass in gram
	zz <- vector(length=0)
	Mass <- vector(length=0)
	MassProtein <-	for(i in 1:length(TotalIntenisty)){
		zz <- (data[,i] * y) / TotalIntenisty[i] * 10^-12
		Mass <- cbind(Mass, zz)
	}
	colnames(Mass) = colnames(data)
	# calculate moles
	Mol <- Mass/(mol*1000)
	Copy <- Mol * 6.022140857*10^23
}

###
# FIT WIBULL
###

##x=timeN;y=P[j,naive]
fit.weibull <- function(x,y) {
    y = as.numeric(y)
    x = as.numeric(x)
    jj <- which( x==0 | y>0)  ## skip zero measurements if not x==0
    x = x[jj]
    y = y[jj]
    cdf.weibull <- function(x, lambda, k) (1 - exp(-(x/lambda)^k))
    fit = nls( y ~ cdf.weibull(x,lambda,k),
              start=list(lambda=50,k=1),
              lower=list(lambda=0.01,k=0.001),
              algorithm="port")
    summary(fit)
    coef(fit)
    xfit = seq(0,48,1)
    lambda0 = coef(fit)["lambda"]
    k0 = coef(fit)["k"]
    yfit = cdf.weibull(xfit, lambda0, k0)
    t50 = lambda0*log(2)**(1/k0)  ## time at 50%
    t50
    list(x=xfit, y=yfit, t50=t50)
}

###Functions for Proteomic Analysis Workflow
library(compiler)
library(zipfR)
library(ggplot2)
library(gplots)
###LOAD DATA
load.data <- function(x, y=getwd()){
  setwd(y)
  z <- read.delim(x, header=T, stringsAsFactor=F)
}

load.data2 <- function(x, y=getwd()){
  setwd(y)
  rows5.df <- read.table(x, header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "\"", comment.char = "#", nrows=3)
  classes <- sapply(rows5.df, class)
  #classes <- ifelse(classes == 'integer', 'character', classes)
  #classes <- ifelse(classes == 'logical', 'character', classes)
  read.table(x, header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "\"", comment.char = "#", colClasses = classes)
}

###SAVE dataframe
safe.data <- function(x, name="data.txt"){
  write.table(x, name, sep="\t", row.names = F)
}

#groups.to.numbers
groupstonumbers <- function (x){
  number.vector <- c(1:length(x))
  i=1
  for(z in unique(x)){
    number.vector[groups==z] <- i
    i=i+1
  }
  return(number.vector)
}

#groups.to.colors
groupstocolors <- function (x){
  color.vector <- c(1:length(x))
  group.color <- rainbow(length(unique(x)))
  i=1
  for(z in unique(x)){
    color.vector[groups==z] <- group.color[i]
    i=i+1
  }
  return(color.vector)
}


###remove only identified by side, contaminants, reverse
filter.reverse.contaminant.only <- function(x, rev=T, con=T, only=T){
  if(rev==T){x <- x[x$Reverse!="+",]}
  if(con==T & sum(colnames(x)=="Potential.contaminant")==1){x <- x[x$Potential.contaminant!="+",]}
  if(con==T & sum(colnames(x)=="Contaminant")==1){x <- x[x$Contaminant!="+",]}
  if(only==T){x <- x[x$Only.identified.by.site!="+",]}
  return(x)
}

###FILTER VALID VALUES
#Note: column names have to be group names
filter.valid.values <- function(x, number.valid=3, groups=colnames(x)){
  #define unique groups
  unique.groups <- unique(groups)
  #create dataframe containing number of valid values for each group
  a <- data.frame(matrix(nrow=nrow(x), ncol=0))
  for(z in unique.groups){
    a <- cbind(a, apply(x[,groups==z]>0, 1, sum, na.rm=T))
  }
  names(a) <- unique.groups
  #check for each row if at least one group fullfills the number of valid values specified
  valid <- apply(a>=number.valid, 1, sum)>=1
  #filter data
  x <- x[valid,]
}

###logarithm
logarithm <- function(x, logbase=2, replace.inf = NA){
  if(logbase==2){x <- log2(x)}
  if(logbase==10){x <- log10(x)}
  x[x==-Inf]  <- replace.inf
  return(x)
}

###IMPUTATION
imputation <- function(x, width=0.3, downshift=1.8, bycolumn=T){
  if(bycolumn==T){
    for(i in 1:ncol(x)){
      x[,i][is.na(x[,i])] <- rnorm(sum(is.na(x[,i])), mean(as.numeric(x[,i]), na.rm=T) - downshift*sd(as.numeric(x[,i]), na.rm=T), width*sd(as.numeric(x[,i]), na.rm=T))
    }
  }else{
    x[is.na(x)] <- rnorm(sum(is.na(x)), mean(as.numeric(x), na.rm=T) - downshift*sd(as.numeric(x), na.rm=T), width*sd(as.numeric(x), na.rm=T))
  }
  return(x)
}

###T-test one vs all
Ttest.one.vs.all <- function(x, groups=colnames(x)){
  #define dataframes
  pvalue.df <- data.frame(matrix(nrow=nrow(x), ncol=0))
  difference.df <- data.frame(matrix(nrow=nrow(x), ncol=0))
  unique.groups <- unique(groups)

  ###test one celltype vs all other cell types
  for(z in unique.groups){
    #define ttest groups
    ttestgroups <- groups==z
    #compute pvalue and difference
    stat.pvalue <- apply(x, 1, function(y){t.test(y[!ttestgroups], y[ttestgroups], var.equal=F)$p.value})
    stat.est.means <- apply(x, 1, function(y){t.test(y[!ttestgroups], y[ttestgroups], var.equal=F)$estimate})
    stat.mean <- apply(stat.est.means, 2, diff)
    print(paste("Celltype", z, sep=":"))
    #write dataframe
    pvalue.df <- cbind(pvalue.df, stat.pvalue)
    difference.df <- cbind(difference.df, stat.mean)
  }
  #renames columns
  colnames(pvalue.df) <- paste("t.test.pvalue", unique.groups, "vs.ALL", sep=".")
  colnames(difference.df) <- paste("t.test.difference", unique.groups, "vs.ALL", sep=".")
  #combine difference and pvalue
  result <- cbind(difference.df, pvalue.df)
  return(result)
}

###T-test combinations
Ttest.combinations <- function(x, groups=colnames(x)){
  #define dataframes
  pvalue.df2 <- data.frame(matrix(nrow=nrow(x), ncol=0))
  difference.df2 <- data.frame(matrix(nrow=nrow(x), ncol=0))
  #define combinations
  unique.groups <- unique(groups)
  celltype.combinations <- combn(unique.groups, 2)
  ###test one celltype vs another cell type
  for(i in 1:ncol(celltype.combinations)){
    #define ttest groups
    ttestgroup1 <- groups==celltype.combinations[1,i]
    ttestgroup2 <- groups==celltype.combinations[2,i]
    #compute pvalue and difference
    stat.pvalue <- apply(x, 1, function(y){t.test(y[ttestgroup1], y[ttestgroup2], var.equal=F)$p.value})
    stat.est.means <- apply(x, 1, function(y){t.test(y[ttestgroup1], y[ttestgroup2], var.equal=F)$estimate})
    stat.mean <- apply(stat.est.means, 2, diff)
    print(paste(i, celltype.combinations[1,i], celltype.combinations[2,i], sep="; "))
    #write dataframe
    pvalue.df2 <- cbind(pvalue.df2, stat.pvalue)
    difference.df2 <- cbind(difference.df2, stat.mean)
  }
  #renames columns
  colnames(pvalue.df2) <- paste("t.test.pvalue", celltype.combinations[1,],"vs", celltype.combinations[2,], sep=".")
  colnames(difference.df2) <- paste("t.test.difference", celltype.combinations[1,],"vs", celltype.combinations[2,], sep=".")
  #combine difference and pvalue
  result <- cbind(difference.df2, pvalue.df2)
  return(result)
}

###VOLCANOPLOT
volcanoplot <- function(x, y, difference.cut=2.5, pvalue.cut=0.01, gene.names=NULL, title=NULL){
  significance <- x^2>difference.cut^2 & y < pvalue.cut
  #sig points
  diff.out <- x[significance]
  pvalue.out <- -log10(y[significance])
  #non sig points
  diff <- x[!significance]
  pvalue <- -log10(y[!significance])
  #gene.names of sig points
  gene.names <- gene.names[significance]

  #scatterplot
  plot(diff, pvalue,
       pch= 19,
       col="#99999920",
       xlab="Fold Change [log2]",
       ylab="p-value [-log10]",
       main=title,
       xlim=c(min(x)-1,max(x)+1),
       ylim=c(0,max(-log10(y))+1)
  )
  #points in sig marked blue
  if(length(diff.out)!=0){
  points(diff.out,pvalue.out,
         pch= 19,
         col="#377EB8")
  }
  #label sig points with gene.name
  if(length(gene.names)!=0){
  text(diff.out,pvalue.out,
       labels=gene.names,
       col="#377EB8", cex=.7, pos=3)
  }
}


###PCA
pca.j <- function(x, mycolor="black", combinations.number=3,
                  sample.names=c(1:ncol(x)), sample.col = "black", gene.names=c(1:nrow(x)),
                  gene.names.filter=c(T, rep(F, nrow(x)-1)), gene.names.col="red"){

  pca.prot <- prcomp(t(x), scale=F)
  variance.percent <- round(pca.prot$sdev^2/sum(pca.prot$sdev^2)*100, 1)

  #how many combinations
  combinations <- combn(c(1:combinations.number),2)

  for(k in 1:ncol(combinations)){
    par(mfcol = c(1, 2))
    ###plot pca components
    i=combinations[1,k]
    j=combinations[2,k]
    x.axis.lab <- paste(colnames(pca.prot$x)[i], " (", variance.percent[i], " %)", sep="")
    y.axis.lab <- paste(colnames(pca.prot$x)[j], " (", variance.percent[j], " %)", sep="")
    plot(pca.prot$x[,i], pca.prot$x[,j],
         col=mycolor, cex=1,
         xlab=x.axis.lab, ylab=y.axis.lab, main="PCA", type="n")
    text(pca.prot$x[,i], pca.prot$x[,j], labels=sample.names, col=sample.col)
    abline(v=0, col="gray", lty = 9)
    abline(h=0, col="gray", lty = 9)
      ###plot pca loadings
    plot(pca.prot$rotation[,i], pca.prot$rotation[,j], pch=21, cex=.5, xlab=paste("PC",i, sep=""), ylab=paste("PC",j, sep=""), main="Loadings", col="#99999980")
    abline(v=0, col="gray", lty = 9)
    abline(h=0, col="gray", lty = 9)
    text(pca.prot$rotation[gene.names.filter,i], pca.prot$rotation[gene.names.filter,j], gene.names[gene.names.filter], col=gene.names.col, cex=1)
  }
}


###Z-SCORE
z.score <- function(x, byrow=T){
  if(byrow==F) x <- t(x)
  z.mean <- apply(x, 1, mean, na.rm=T)
  z.sd <- apply(x, 1, sd, na.rm=T)
  x <- (x-z.mean)/z.sd
  if(byrow==F) x <- t(x)
  return(x)
}

###HEATMAP
heatmap.j <- function(x, gradient=T, range.white=3, ColSideColors="black", colv=T, rowv=T){

  #library("gplots")

  if(gradient==T){
  my_palette <- colorRampPalette(c("gray", "white", "yellow", "red"))(n = 399)
  a <- max(x, na.rm=T)
  b <- min(x[x!=0], na.rm=T)
  c <- median(x[x!=0], na.rm=T)

  col_breaks = c(seq(0,1,length=100),  # for gray
                 seq(1,b,length=100),  # for white
                 seq(b,c,length=100),              # for yellow
                 seq(c,a,length=100)) # for red
  }
  if(gradient==F){
    my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
    a <- max(x, na.rm=T)
    b <- min(x[x!=0], na.rm=T)

    col_breaks = c(seq(b,-range.white,length=100),  # for blue
                   seq(-range.white,range.white,length=100),              # for white
                   seq(range.white,a,length=100))  # for red
  }
  heatmap.2(x,
            col = my_palette,
            breaks=col_breaks,
            margins = c(5, 5),
            trace = "none",
            xlab = "Cell types",
            ylab = "Gene Names",
            #lmat=rbind(c(4,3,3), c(2,1,1), c(2,1,1)),
            #lhei=c(1,2,3),
            #lwid=c(1,2,3),
            scale = c("none"),
            #symbreaks = min(18, na.rm=TRUE),
            na.color="grey",
            cexRow = .9, cexCol = .9,
            main = "T Cell Clones",
            dendrogram = "both",
            density.info = "none",
            Colv = colv,
            Rowv = rowv,
            ColSideColors=ColSideColors,
            # block sepration
            #colsep = c(1:ncol(x)),
            #rowsep = c(1:nrow(x)),
            #sepcolor="white",
            #sepwidth=c(0.01,0.01),
            #cellnote=round(x, 1),
            #notecol="black",
            #notecex=.7
  )
}

###t test s0
slow_s0.ttest <- function(x, y, s0, alternative="two.sided", na.rm=F)
{
  xval <- as.numeric(x)
  # remove NaN
  xlen <- length(xval)
  xmean <- mean(xval)

  yval <- as.numeric(y)
  # remove NaN
  ylen <- length(yval)
  ymean <- mean(yval)

  s  <- sqrt((sum((xval-xmean)^2) + sum((yval-ymean)^2)) * (1/xlen+1/ylen) / (xlen+ylen-2))
  df <- xlen + ylen - 2
  diff <- xmean - ymean
  df <- xlen + ylen - 2

  t <- diff/s
  pval <- 0.5 * Rbeta(df/(df + t^2), df/2, 0.5)
  t.s0 <- diff/(s + s0)
  pval.s0 <- 0.5 * Rbeta(df/(df + t.s0^2), df/2, 0.5)

  if (alternative == "less") { # go to the correct alternative (currently greater)
    pval <- 1 - pval
    pval.s0 <- 1 - pval.s0
  } else if (alternative == "two.sided") {
    pval <- 2 * min(pval, 1 - pval)
    pval.s0 <- 2 * min(pval.s0, 1 - pval.s0)
  }

  c(
    "df"= df,
    "t" = t,
    "p.value" = pval,
    "t.s0" = t.s0,
    "p.value.s0" = pval.s0
  )
}
s0.ttest <- cmpfun(slow_s0.ttest)

###fdr controlled s0
slow_fdr.controlled.ttest <- function(df, col1, col2, s0=1, permutations=125)
{
  # build up the result data frame
  all.cols <- c(col1,col2)
  ncols <- length(all.cols)
  d.experiments <- as.matrix(df[,c(col1,col2)])
  p.values.for.fdr <- cbind(d.experiments, data.frame(
    Index                    = 1:nrow(df),
    positive                = 1,
    negative              = 0
  ))

  # build up the test-matrix
  for (permutation in 1:permutations)
  {
    rnd.col <- rep(NaN, ncols)
    i.col1 <- 1
    rnd.col1 <- sample(col1)
    i.col2 <- 1
    rnd.col2 <- sample(col2)
    for (i in 1:ncols) {
      if ((i%%2 == 0 && i.col1 <= length(col1)) || i.col2 > length(col2)) {
        rnd.col[i] <- rnd.col1[i.col1]
        i.col1 <- i.col1 + 1
      } else if ((i%%2 == 1 && i.col2 <= length(col2)) || i.col1 > length(col1)) {
        rnd.col[i] <- rnd.col2[i.col2]
        i.col2 <- i.col2 + 1
      }
    }

    d.random <- df[,rnd.col]
    colnames(d.random) <- c(col1,col2)
    p.values.for.fdr <- rbind(p.values.for.fdr, cbind(d.random, data.frame(
      Index                    = rep(-1, nrow(d.random)),
      positive                = 0,
      negative              = 1
    )))
  }

  # calculate the statistic
  p.values.for.fdr <- cbind(
    p.values.for.fdr,
    t(apply(p.values.for.fdr, 1, function(x) s0.ttest(x[col1], x[col2], s0=s0)))
  )

  # calculate q-value on s0
  p.values.for.fdr <- p.values.for.fdr[order(p.values.for.fdr$p.value.s0),]
  p.values.for.fdr$positive <- cumsum(p.values.for.fdr$positive)
  p.values.for.fdr$negative <- cumsum(p.values.for.fdr$negative)
  p.values.for.fdr$q.value <- apply(p.values.for.fdr, 1, function(x) min(1, x["negative"]/x["positive"]/permutations))
  res <- p.values.for.fdr[which(p.values.for.fdr$Index != -1),]

  # smooth out the q-values - i.e. reset the q-value to that of the worst p-value with a better q-value
  for (i in 1:length(res$q.value))
  {
    qval <- res$q.value[i]
    res$q.value[which(res$q.value[1:i] > qval)] <- qval
  }

  # done, cleanup memory and return the result
  res <- res[order(res$Index),]
  res <- res[,!(names(res) %in% c("Index", "positive", "negative", col1, col2))]
  return(res)
}
fdr.controlled.ttest <- cmpfun(slow_fdr.controlled.ttest)

#interaction by genename
interaction.stringDB.gene.names <- function(gene.names=gene.names, confidence.level=400){
  #reomve doubles IDs
  gene.names <- gene.names[!duplicated(gene.names)]
  #reomve ""
  gene.names <- gene.names[!gene.names==""]
  #reomve "2-sep"
  gene.names <- gene.names[!grepl("^[0-9]-Sep$", gene.names)]
  #sort gene.names
  gene.names <- sort(gene.names)
  #safe gene.names in different variable
  gene.names.original <- gene.names
  #remove double genenames ("...;...")
  gene.names <- sub(";.*$", "", gene.names)
  gene.names <- sub("-[0-9]{1,2}$", "", gene.names)
  #get stringIDs
  stringIDs = string_db$mp(gene.names)
  #get interaction data
  protein.interactions <- string_db$get_interactions( stringIDs )
  #add gene names to interaction data
  for(i in 1:length(gene.names)){
    protein.interactions$from.gene.names[protein.interactions$from == stringIDs[i]] <- gene.names.original[i]
    protein.interactions$to.gene.names[protein.interactions$to == stringIDs[i]] <- gene.names.original[i]
  }
  #filter protein interaction data
  protein.interactions.filtered <- protein.interactions[protein.interactions$experimental>=confidence.level,]
  #result
  return(protein.interactions.filtered[12:3])
}


#input needs to be a dataframe containing four columns (id, id.ori, gene, gene.ori)
interaction.stringDB.proteinIDs <- function(proteinIDs.df, confidence.level=400, map.col = "id"){

  proteinIDs.df$id <- sub(";.*$", "", proteinIDs.df$id)
  proteinIDs.df$id <- sub("-[0-9]$", "", proteinIDs.df$id)
  proteinIDs.df$gene <- sub(";.*$", "", proteinIDs.df$gene)
  head(proteinIDs.df)
  #get stringIDs
  stringIDs.df = string_db$map(proteinIDs.df, map.col, removeUnmappedRows = F)
  not.mapped <- stringIDs.df[is.na(stringIDs.df$STRING_id),]

  stringIDs <- stringIDs.df$STRING_id[!duplicated(stringIDs.df$STRING_id) & !is.na(stringIDs.df$STRING_id)]
  gene.names.original <- stringIDs.df$gene.ori[!duplicated(stringIDs.df$STRING_id) & !is.na(stringIDs.df$STRING_id)]
  #get interaction data
  protein.interactions <- string_db$get_interactions( stringIDs )
  #add gene names to interaction data
  for(i in 1:length(stringIDs)){
    protein.interactions$from.gene.names[protein.interactions$from == stringIDs[i]] <- gene.names.original[i]
    protein.interactions$to.gene.names[protein.interactions$to == stringIDs[i]] <- gene.names.original[i]
  }
  #filter protein interaction data
  protein.interactions.filtered <- protein.interactions[protein.interactions$experimental>=confidence.level,]
  #result
  return(list(A=protein.interactions.filtered[12:3], B=not.mapped))
}

#plot interaction network
plot.network <- function(network.data, vertex.thickness=.5, edge.thickness=1, layout=layout.fruchterman.reingold){
  network.x <- graph.data.frame(network.data, directed=F)
  V(network.x)$size <- degree(network.x)*vertex.thickness
  E(network.x)$width <- scale(network.data[,3]) + abs(min(scale(network.data[,3])))*edge.thickness

  par(mai=c(0,0,1,0))
  plot(network.x,
       layout=layout,
       main="Network",
       vertex.label.dist=0.5,
       vertex.frame.color="blue",
       vertex.label.color="black",
       vertex.label.font=2,
       vertex.label=V(network.x)$name,
       vertex.label.cex=1
  )
}

