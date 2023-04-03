##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##


# This file is supposed to run from the root Playground folder
if(basename(getwd()) != "omicsplayground") {
  stop("This file is supposed to run from the root Playground folder")
}

options(Ncpus=8L)
options(repos = c(REPO_NAME = "https://cloud.r-project.org/"))

if(1) {
    ## Speed up installation using binary packages from RStudio. Works only for 20.04 LTS !!!
    options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])))
    source("https://docs.rstudio.com/rspm/admin/check-user-agent.R")
    options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/__linux__/jammy/latest"))
}

install.packages("devtools")
install.packages("BiocManager")
install.packages("renv")
LOCAL.PKGS <- sub("_.*","",dir("extdata/packages"))
## LOCAL.PKGS
INSTALLED.PKGS = c("devtools","BiocManager","renv")

install.pkg <- function(pkg, force=FALSE) {
    if(force || !pkg %in% rownames(installed.packages())) {
        if(pkg %in% LOCAL.PKGS) {
            ## if available locally, we install local version
            cat("installing",pkg,"from local folder...\n")
            pkg1 = dir("extdata/packages",pattern=paste0(pkg,"_"),full.names=TRUE)
            try(install.packages(pkg1,repos=NULL,type="source",force=force))
        } else {
            cat("installing",pkg,"from CRAN/BioConductor...\n")
            try(BiocManager::install(
                 pkg, dependencies=NA,force=force,ask=FALSE,update=FALSE))
            if(!require(pkg, character.only=TRUE)) {
                cat("retrying to install",pkg,"from CRAN...\n")
                try(install.packages(
                    pkg, dependencies=NA, force=force, ask=FALSE, update=FALSE))
            }
        }
    } else {
        cat("package",pkg,"already installed\n")
    }
    if(require(pkg, character.only=TRUE)) {
        INSTALLED.PKGS <<- unique(c(INSTALLED.PKGS, pkg))
    } 
}
install.pkgs <- function(pkgs, force=FALSE) {
    pkgs <- sort(unique(pkgs))
    for(pkg in pkgs) install.pkg(pkg, force=force)
}
remove.pkg <- function(pkg) {
    if(pkg %in% rownames(installed.packages())) try(remove.packages(pkg))
    INSTALLED.PKGS <- setdiff(INSTALLED.PKGS, pkg)
}
remove.pkgs <- function(pkgs, force=FALSE) {
    pkgs <- sort(unique(pkgs))
    for(pkg in pkgs) {
        cat("removing",pkg,"\n")
        remove.pkg(pkg)
        INSTALLED.PKGS <<- setdiff(INSTALLED.PKGS, pkg)
    }
}
install.github <- function(repo, force=FALSE) {
    pkg <- sub(".*/","",repo)
    if(!require(pkg, character.only=TRUE) || force) {
      devtools::install_github(repo, upgrad="never", build_vignettes=FALSE, force=TRUE)
    } else {
        cat("package",repo,"already installed\n")      
    }
    if(require(pkg, character.only=TRUE)) {
        INSTALLED.PKGS <<- unique(c(INSTALLED.PKGS, pkg))
    } 
}
autoscan.pkgs <- function() {
    rfiles1 <- system("find components -name \\*.r", intern=TRUE)
    rfiles2 <- system("find components -name \\*.R", intern=TRUE)
    rfiles <- paste(c(rfiles1,rfiles2),collapse=" ")
    pkg1 <- system(paste("grep '::' ",rfiles), intern=TRUE)
    pkg2 <- system(paste("grep 'require(' ",rfiles), intern=TRUE)
    pkg3 <- system(paste("grep 'library(' ",rfiles), intern=TRUE)
    pkg <- c(pkg1,pkg2,pkg3)
    pkg <- grep("message|dbg|cat",pkg,value=TRUE,invert=TRUE)
    pkg <- gsub(".*[rR]:","",pkg)  ## strip filename
    pkg <- grep("^[#]",pkg,invert=TRUE,value=TRUE)  ## no comments
    pkg <- trimws(pkg)
    pkg <- gsub("[:\"]","",gsub(".*[ ,\\(\\[]","",gsub("::.*","::",pkg)))
    pkg <- gsub("\\).*","",gsub(".*require\\(","",pkg))
    pkg <- gsub("\\).*","",gsub(".*library\\(","",pkg))
    pkg <- grep("[=#/*'\\]",pkg,value=TRUE,invert=TRUE)  ## skip commented out
    pkg <- grep("^[a-z]",pkg,value=TRUE,ignore.case=TRUE)  ## skip commented out
    pkg <- grep("[a-z.]*",pkg,value=TRUE,ignore.case=TRUE)  ## skip commented out
    pkg <- setdiff(pkg,c(""))
    pkg <- sort(unique(pkg))
    pkg
}

##---------------------------------------------------------------------
## Install basic shiny packages
##---------------------------------------------------------------------

base.pkg = c("shiny","flexdashboard","shinydashboard", "shinyBS","systemfonts","shinyjs",
             "shinydashboardPlus",'R.utils','shinythemes',"shinybusy","shinycssloaders",
             "shinyWidgets")
install.pkgs(base.pkg, force=FALSE)

##---------------------------------------------------------------------
## Automatically scan all used packages and install
##---------------------------------------------------------------------

##pkg.used <- autoscan.pkgs()

## We use renv to detect dependencies. Renv is looking for library and
## require statements in the r/R source files.
renv.out <- renv::dependencies(path="components")
#head(renv.out)
pkg.used <- unique(renv.out$Package)

## These packages were not detected by renv::dependencies(). We should
## check if they are actually used or needed.
pkg.extra <- c(
  "BioBase","SingleCellExperiment","preprocessCore",
  "liger","monocle3","bsutils","reshape","waiter","sever",
  "RSpectra","SmartSVA","SILGGM","flashClust","sortable",
  "TCGAbiolinks","TCGAutils","GEOmetadb","Rtsne",
  'wordcloud2', 'optparse', 'docopt',"DT","plotly",
  'kableExtra', 'shinythemes', 'rworldmap',
  "HiddenMarkov","coin","rjags","argparse",
  "RcppParallel", "KEGGgraph", 
  'TxDb.Hsapiens.UCSC.hg19.knownGene',
  'TxDb.Mmusculus.UCSC.mm10.knownGene',
  'listviewer'
)

pkg.used <- c(pkg.used, pkg.extra)
pkg.used <- sort(unique(pkg.used))

## We install these later because some want dependencies alreay
## installed.
pkg.later <- c(
    "gputools","Seurat","EPIC","NNLM","iTALK",
    "fpc","grid","gridGraphics","Rgraphviz", ## "rWordCloud",
    "FastGGM","monocle3","proteus",
    "infercnv","pathview",
    "mygene","diptest","edgeR","DESeq2"
  )

install.pkgs( setdiff(pkg.used,pkg.later) )

length(INSTALLED.PKGS)

##---------------------------------------------------------------------
## reinstall problematics ones
##---------------------------------------------------------------------

##install.pkg("grid", force=TRUE)
install.pkgs(c("gridGraphics","Rgraphviz","fastcluster", "mygene",
               "diptest", "fpc"))
## install webshot and phantomjs (do we need it?)
#install.pkg("webshot")
#webshot::install_phantomjs(force=TRUE)  ## cp to /usr/local/bin !!
#file.copy("~/bin/phantomjs","/usr/local/bin") ## need sudo!!
#file.copy("/root/bin/phantomjs","/usr/local/bin") ## inside docker we are root

devtools::install_version("mnormt", repos="http://cran.us.r-project.org") ## for metap
install.pkgs(c('umap','corrplot','wordcloud','metap','brew'))
install.pkgs(c('monocle','Seurat'))
install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz")
install.packages('https://www.bioconductor.org/packages/3.11/data/annotation/src/contrib/KEGG.db_3.2.4.tar.gz')

##---------------------------------------------------------------------
## Install latest from GITHUB (overwriting any other version)
##---------------------------------------------------------------------
install.github("GfellerLab/EPIC")
##install.github("IOR-Bioinformatics/PCSF", dependencies=TRUE, type="source")
install.github('linxihui/NNLM')
install.github("Coolgenome/iTALK")
install.github('wt2015-github/FastGGM', force=TRUE)
install.github("JohnCoene/waiter")
install.github('JohnCoene/firebase@omics', force=TRUE)

##---------------------------------------------------------------------
## ONLY DEV.MODE (single-cell trajectories)
##---------------------------------------------------------------------
if(1) {
    ## ---- monocle3 (only DEV!!! many install problems in R 3.5.2!!!)
    install.pkgs(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                   'limma', 'S4Vectors', 'SingleCellExperiment',
                   'SummarizedExperiment', 'batchelor'))
    install.github('cole-trapnell-lab/leidenbase')
    install.github('cole-trapnell-lab/monocle3')
}

## proteus
devtools::install_github("bartongroup/Proteus", build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=FALSE)
INSTALLED.PKGS <- c(INSTALLED.PKGS, "Proteus")

##---------------------------------------------------------------------
## make sure LOCAL ones are preferred and overwriting previously
## installed versions
## ---------------------------------------------------------------------
## force=TRUE
if(getUsername.System()=="root") system("sudo apt -y install jags")
install.pkg("rjags", force=TRUE)
install.pkgs(LOCAL.PKGS, force=TRUE)

##---------------------------------------------------------------------
## Install Kaleido for plotly
##---------------------------------------------------------------------

install.pkg('reticulate')
reticulate::install_miniconda()
reticulate::conda_install('r-reticulate', 'python-kaleido')
reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
reticulate::use_miniconda('r-reticulate')

##---------------------------------------------------------------------
## remove unneccessary big packages??
##---------------------------------------------------------------------
BIG.NOTUSED <- c(
    "reactome.db", ## >2GB!!!
    "BH","PCSF","terra",
    "DeMixT", ## purify
    "RNAseqData.HNRNPC.bam.chr14",
    "RSpectra",
    ##"org.Mm.eg.db",
    "tximportData"
    ##"EnsDb.Hsapiens.v86",
    ##"EnsDb.Mmusculus.v79",
    ##"TxDb.Hsapiens.UCSC.hg19.knownGene",  ## need for import
    ##"TxDb.Mmusculus.UCSC.mm10.knownGene"  ## need for import
)
remove.pkgs(BIG.NOTUSED)

if(1) {
    ## Write license file of the used/installed packages
    lisc <- installed.packages(fields = "License")
    sel <- which(lisc[,"Package"] %in% INSTALLED.PKGS)
    lisc1 <- lisc[sel,]
    lisc1 <- lisc1[order(lisc1[,"Package"]),]
    lisc2 <- lisc1[,c("Package","Version","License")]
    ##write.table(lisc2, "RPackageLicenses.txt",sep='\t', quote=FALSE, row.names=FALSE)
    fixstr <- function(s,n=30) {substring(paste0(s,paste(rep(" ",n),collapse='')),1,n) }
    lisc2 <- rbind( c("PACKAGE","VERSION","LICENSE"), lisc2)
    lisc2 <- cbind(fixstr(lisc2[,1],36),fixstr(lisc2[,2],15),fixstr(lisc2[,3],30))
    lisc.text <- apply(lisc2, 1, function(s) paste0(s,collapse=''))
    write(lisc.text, "RPackageLicenses.txt")
##  lisc2[grep("LGPL|AGPL|GPL-3",lisc2[,"License"]),]
  
  
}
