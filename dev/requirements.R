##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##
## NOTE: This file is supposed to run in the folder .../R/
##

options(repos = c(REPO_NAME = "https://cloud.r-project.org/"))

if(0) {
    ## Speed up installation using binary packages from RStudio. Works only for 20.04 LTS !!!
    options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])))
    source("https://docs.rstudio.com/rspm/admin/check-user-agent.R")
    options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/__linux__/focal/latest")) 
}

options(Ncpus=8L)

install.packages("devtools")
install.packages("BiocManager")

LOCAL.PKGS <- sub("_.*","",dir("../ext/packages"))
LOCAL.PKGS

install.pkg <- function(pkg, force=FALSE) {
    if(force || !pkg %in% rownames(installed.packages())) {
        if(pkg %in% LOCAL.PKGS) {
            ## if available locally, we install local version
            cat("installing",pkg,"from local folder...\n")
            pkg1 = dir("../ext/packages",pattern=paste0(pkg,"_"),full.names=TRUE)
            try(install.packages(pkg1,repos=NULL,type="source",force=force))
        } else {
            cat("installing",pkg,"from CRAN/BioConductor...\n")
            try(BiocManager::install(
                 pkg, dependencies=NA,force=force,ask=FALSE,update=FALSE))
            if(!require(pkg, character.only=TRUE)){
                cat("retrying to install",pkg,"from CRAN...\n")
                try(install.packages(
                    pkg, dependencies=NA, force=force, ask=FALSE, update=FALSE))
            }
        }
    } else {
        cat("package",pkg,"already installed\n")
    }
}
install.pkgs <- function(pkgs, force=FALSE) {
    pkgs <- sort(unique(pkgs))
    for(pkg in pkgs) install.pkg(pkg, force=force)
}
remove.pkg <- function(pkg) {
    if(pkg %in% rownames(installed.packages())) try(remove.packages(pkg))
}
remove.pkgs <- function(pkgs, force=FALSE) {
    pkgs <- sort(unique(pkgs))
    for(pkg in pkgs) {
        cat("removing",pkg,"\n")        
        remove.pkg(pkg)
    }
}

autoscan.pkgs <- function() {

    pkg1 <- system("grep '::' *.r *.R ../shiny/boards/*R ../shiny/modules/*R", intern=TRUE)
    pkg2 <- system("grep 'require(' *.r *.R ../shiny/boards/*R ../shiny/modules/*R", intern=TRUE)
    pkg3 <- system("grep 'library(' *.r *.R ../shiny/boards/*R ../shiny/modules/*R", intern=TRUE)    

    pkg <- c(pkg1,pkg2,pkg3)
    pkg <- grep("message|dbg|cat",pkg,value=TRUE,invert=TRUE)

    pkg <- gsub(".*[rR]:","",pkg)  ## strip filename
    pkg <- grep("^[#]",pkg,invert=TRUE,value=TRUE)  ## no comments
    
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

PKG.MANUAL <- c(
    "gputools","Seurat","EPIC","PCSF","NNLM","iTALK",
    "fpc","grid","gridGraphics","Rgraphviz", ## "rWordCloud",
    "shinyparticles","FastGGM","monocle3","proteus",
    ## "infercnv","pathview",
    "mygene","diptest","edgeR","DESeq2")

##---------------------------------------------------------------------
## Install base packages
##---------------------------------------------------------------------

base.pkg = c("shiny","flexdashboard","shinydashboard",
             "shinydashboardPlus",'R.utils','shinythemes')
install.pkgs(base.pkg, force=TRUE)

##---------------------------------------------------------------------
## Automatically scan all used packages and install
##---------------------------------------------------------------------

## pkg.used <- system("grep 'library(\\|require(' *R *r ../shiny/*R ../shiny/modules/*R", intern=TRUE)
## pkg.used <- gsub(".*require\\(|.*library\\(","",pkg.used)
## pkg.used <- gsub("\"|\'|\\).*","",pkg.used)
## pkg.used <- grep("[ ]|quietly",pkg.used,value=TRUE,invert=TRUE)

pkg.used <- autoscan.pkgs() 

pkg.needed <- c('umap', 'corrplot', 'wordcloud', 'wordcloud2', 'optparse', 'docopt',
                'kableExtra', 'randomForest', 'rhdf5', 'qgraph', 'psych',
                'ggVennDiagram', 'shinythemes', 'shinybusy', 'beepr',
                'rworldmap', 'sever', 'WGCNA', 'DGCA', 'ggforce')

pkg.used <- c(pkg.used, pkg.needed)
pkg.used <- sort(unique(pkg.used))
install.pkgs( setdiff(pkg.used,PKG.MANUAL) )

r.pkg <- c('TxDb.Hsapiens.UCSC.hg19.knownGene',
           'TxDb.Mmusculus.UCSC.mm10.knownGene')
install.pkgs(r.pkg)

##---------------------------------------------------------------------
## reinstall problematics ones
##---------------------------------------------------------------------

##install.pkg("grid", force=TRUE)
install.pkgs(c("gridGraphics","Rgraphviz","fastcluster", "mygene",
               "diptest", "fpc", "webshot"))
webshot::install_phantomjs(force=TRUE)  ## cp to /usr/local/bin !!
file.copy("~/bin/phantomjs","/usr/local/bin") ## need sudo!!

devtools::install_version("mnormt", version="1.5-7", repos="http://cran.us.r-project.org")
install.pkgs(c('umap','corrplot','wordcloud','metap','brew'))
install.pkgs(c('monocle','Seurat'))
install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz")

##---------------------------------------------------------------------
## Install latest from GITHUB
##---------------------------------------------------------------------
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
##devtools::install_github("IOR-Bioinformatics/PCSF", dependencies=TRUE, type="source")
devtools::install_github('linxihui/NNLM')
devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE, force=TRUE)
## devtools::install_github('adymimos/rWordCloud', force=TRUE)
## remotes::install_github("dreamRs/shinyparticles")
remotes::install_github("dreamRs/particlesjs")
remotes::install_github("JohnCoene/waiter")

##---------------------------------------------------------------------
## ONLY DEV.MODE (single-cell trajectories)
##---------------------------------------------------------------------
if(1) {
    ## ---- monocle3 (only DEV!!! many install problems in R 3.5.2!!!)
    install.pkgs(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                   'limma', 'S4Vectors', 'SingleCellExperiment',
                   'SummarizedExperiment', 'batchelor'))
    devtools::install_github('cole-trapnell-lab/leidenbase')
    devtools::install_github('cole-trapnell-lab/monocle3')
}

## proteus
devtools::install_github("bartongroup/Proteus", build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=FALSE)

##---------------------------------------------------------------------
## Seurat needs to be downgraded and dependencies to installed...
##---------------------------------------------------------------------
## install.pkgs(c("ROCR", "mixtools", "lars", "ica", "tsne", "ape", "dtw", "SDMTools", "ggridges", "fitdistrplus", "doSNOW","diffusionMap","fpc","hdf5r"))
## install.pkgs(c('cowplot', 'Rtsne', 'pbapply', 'RANN', 'dplyr', 'irlba', 'plotly', 'Hmisc', 'tidyr', 'metap', 'lmtest', 'png', 'reticulate', 'RcppEigen', 'RcppProgress'))
## install.packages("../ext/packages/Seurat_v2.3.3.tar.gz",repos=NULL,type="source")  ## old version

##---------------------------------------------------------------------
## make sure local ones are preferred
##---------------------------------------------------------------------
##force=TRUE

install.pkgs(c("RcppParallel"))
devtools::install_github('wt2015-github/FastGGM', force=TRUE)
install.pkgs(c("HiddenMarkov","coin","rjags","future","argparse"))
install.pkg("rjags", force=TRUE)
##install.packages("../ext/packages/infercnv_1.1.3mod.tar.gz",repos=NULL,type="source")  ## old version
BiocManager::install("infercnv")
install.pkgs(c("KEGGREST","pathview"), force=TRUE)

##---------------------------------------------------------------------
## Install Kaleido for plotly
##---------------------------------------------------------------------

install.packages('reticulate')
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


if(0) {
    
    pkg1 <- system("grep '::' *.r *.R ../shiny/boards/*R ../shiny/modules/*R", intern=TRUE)
    pkg2 <- system("grep 'require(' *.r *.R ../shiny/boards/*R ../shiny/modules/*R", intern=TRUE)
    pkg3 <- system("grep 'library(' *.r *.R ../shiny/boards/*R ../shiny/modules/*R", intern=TRUE)    

    pkg <- c(pkg1,pkg2,pkg3)
    pkg <- grep("message|dbg|cat",pkg,value=TRUE,invert=TRUE)
    
    pkg <- gsub("[:\"]","",gsub(".*[ ,\\(\\[]","",gsub("::.*","::",pkg)))
    pkg <- gsub("\\).*","",gsub(".*require\\(","",pkg))
    pkg <- gsub("\\).*","",gsub(".*library\\(","",pkg))    
    pkg <- grep("[=#/*'\\]",pkg,value=TRUE,invert=TRUE)

    pkg <- unique(pkg)
    
    lisc <- installed.packages(fields = "License")
    sel <- which(lisc[,"Package"] %in% pkg)
    ##pkg2 <- c(lisc[sel,"Package"], lisc[sel,"Imports"], lisc[sel,"LinkingTo"])
    pkg2 <- c(lisc[sel,"Package"])
    pkg2 <- setdiff(pkg2, c("",NA))
    pkg2 <- unique(sub("[,]","",unlist(lapply(pkg2, function(p) strsplit(p, split='[ ,]')))))

    pkg2 <- intersect(pkg2,lisc[,"Package"])
    lisc1 <- lisc[which(lisc[,"Package"] %in% pkg2),]
    lisc1 <- lisc1[,c(1,3,10)]
    lisc1 <- lisc1[order(lisc1[,"Package"]),]
    write.table(lisc1, "RPackageLicenses.txt",sep='\t', quote=FALSE, row.names=FALSE)

    fixstr <- function(s,n=30) {substring(paste0(s,paste(rep(" ",n),collapse='')),1,n) }
    lisc2 <- cbind(fixstr(lisc1[,1],36),fixstr(lisc1[,2],15),fixstr(lisc1[,3],30))
    apply(lisc2, 1, function(s) cat(" -",paste0(s),'\n'))
        
    lisc1[grep("AGPL|GPL-3",lisc1[,"License"]),]

    ## Using binary packages from RStudio
    options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])))
    source("https://docs.rstudio.com/rspm/admin/check-user-agent.R")
    options(repos = c(REPO_NAME = "https://cloud.r-project.org/"))
    options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/latest"))
    options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/__linux__/focal/latest"))
    options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/__linux__/jammy/latest"))
    ##options(repos = c(REPO_NAME = "https://stat.ethz.ch/CRAN/"))
    ##options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/latest"))
    remove.packages("shiny")
    install.packages("shiny")
    packageVersion('shiny')
    
    remove.packages("dplyr")
    install.packages("dplyr")
    packageVersion('dplyr')

    remove.packages("stringi")
    install.packages("stringi")
    packageVersion('stringi')

    BiocManager::install("dplyr", dependencies=NA, ask=FALSE, update=FALSE, force=TRUE)
    install.pkg("dplyr")
    
}
