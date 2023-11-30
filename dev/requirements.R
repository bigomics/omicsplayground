##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


# This file is supposed to run from the root Playground folder
if(basename(getwd()) != "omicsplayground") {
  stop("This file is supposed to run from the root Playground folder")
}

options(Ncpus = 8L)
options(timeout = 99999)  ## download time.out
options(repos = c(REPO_NAME = "https://cloud.r-project.org/"))

## Speed up installation using binary packages from RStudio.
if(grepl("linux",R.version["os"])) {
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
      devtools::install_github(repo, upgrade="never", build_vignettes=FALSE, force=force)
    } else {
        cat("package",repo,"already installed\n")
    }
    if(require(pkg, character.only=TRUE)) {
        INSTALLED.PKGS <<- unique(c(INSTALLED.PKGS, pkg))
    }
}

##---------------------------------------------------------------------
## Install basic shiny packages
##---------------------------------------------------------------------

base.pkg = c("shiny","flexdashboard","shinydashboard", "shinyBS","systemfonts","shinyjs",
             "shinydashboardPlus",'R.utils','shinythemes',"shinybusy","shinycssloaders",
             "shinyWidgets")

##---------------------------------------------------------------------
## Automatically scan all used packages and install
##---------------------------------------------------------------------

## We use renv to detect dependencies. Renv is looking for library and
## require statements in the r/R source files.
##install.packages("renv")
cat("RENV:: building dependencies...\n")
renv.out <- renv::dependencies(path="components", root=getwd(), errors="ignored")
pkg.used <- unique(renv.out$Package)
cat("RENV:: done!\n")

## These packages were not detected by renv::dependencies(). We should
## check if they are actually used or needed.
pkg.extra <- c(
  "BioBase","SingleCellExperiment","preprocessCore",
  "liger","monocle3","bsutils","reshape","waiter","sever",
  "RSpectra","SmartSVA","SILGGM","flashClust","ggrepel", "ComplexHeatmap",
  "TCGAbiolinks","TCGAutils","GEOmetadb","Rtsne", "seriation","sortable",
  'wordcloud2', 'optparse', 'docopt',"DT","plotly",
  'kableExtra', 'shinythemes', 'rworldmap',"e1071","mixOmics",
  "HiddenMarkov","coin","rjags","argparse",
  "RcppParallel", "KEGGgraph", "svgPanZoom",
  'TxDb.Hsapiens.UCSC.hg19.knownGene',
  'TxDb.Mmusculus.UCSC.mm10.knownGene',
  'listviewer','SBGNview','org.Hs.eg.db','DeMixT',
  'svgPanZoom','rhdf5','monocle','mygene',
  'iheatmapr','RcppZiggurat','Rfast','BH','topGO', 'survcomp','plsRcox',
  'blastula','shinytest2'
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


##---------------------------------------------------------------------
## start install
##---------------------------------------------------------------------

install.pkgs(base.pkg, force=FALSE)
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
install.github('linxihui/NNLM')
install.github("Coolgenome/iTALK")
install.github('wt2015-github/FastGGM', force=TRUE)
install.github("JohnCoene/waiter")
install.github('JohnCoene/firebase@omics', force=TRUE)
install.github('JohnCoene/bsutils')
install.github('bigomics/bigdash', force=TRUE)
install.github('bigomics/bigLoaders')
install.github('bigomics/PCSF', force=TRUE)
install.github('bigomics/shinyChatR')
install.github('bigomics/fgsea')
install.github('ropensci/iheatmapr')
install.github('rstudio/bslib@v0.5.1',dependencies=FALSE)
install.github('rstudio/htmltools',dependencies=FALSE)
install_github('bigomics/biomaRt',dependencies=FALSE)

##---------------------------------------------------------------------
## ONLY DEV.MODE (single-cell trajectories)
##---------------------------------------------------------------------
if(1) {
    ## ---- monocle3 (only DEV!!! many install problems in R 3.5.2!!!)
    install.pkgs(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                   'limma', 'S4Vectors', 'SingleCellExperiment',
                   'SummarizedExperiment', 'batchelor'))
    install.github('cole-trapnell-lab/leidenbase')
    #install.github('cole-trapnell-lab/monocle3')
}

## proteus
devtools::install_github("bartongroup/Proteus", build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=FALSE)

#install rliger and dependencies (riverplot)
devtools::install_github("cran/riverplot")
BiocManager::install("rliger")

# install maptools
devtools::install_github("cran/maptools")

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
reticulate::install_miniconda(force=TRUE)
reticulate::conda_install('r-reticulate', 'python-kaleido')
reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
reticulate::use_miniconda('r-reticulate')

##---------------------------------------------------------------------
## remove unneccessary big packages??
##---------------------------------------------------------------------
BIG.NOTUSED <- c(
    "reactome.db", ## >2GB!!!
    "terra",
    ## "DeMixT", ## purify
    "RNAseqData.HNRNPC.bam.chr14",
    ## "RSpectra",  ## ???
    ##"org.Mm.eg.db",
    "tximportData"
    ##"EnsDb.Hsapiens.v86",
    ##"EnsDb.Mmusculus.v79",
    ##"TxDb.Hsapiens.UCSC.hg19.knownGene",  ## need for import
    ##"TxDb.Mmusculus.UCSC.mm10.knownGene"  ## need for import
)
remove.pkgs(BIG.NOTUSED)

## --------------------------------------------------
## Write license file of the used/installed packages
## --------------------------------------------------
if(1) {
    ## Write license file of the used/installed packages
    lisc <- installed.packages(fields = "License")
    sel <- which(lisc[,"Package"] %in% INSTALLED.PKGS)
    lisc1 <- lisc
    lisc1 <- lisc[sel,]    
    lisc1 <- lisc1[order(lisc1[,"Package"]),]
    lisc1 <- lisc1[!duplicated(lisc1[,"Package"]),]
    lisc2 <- lisc1[,c("Package","Version","License")]
    ##write.table(lisc2, "RPackageLicenses.txt",sep='\t', quote=FALSE, row.names=FALSE)
    fixstr <- function(s,n=30) {substring(paste0(s,paste(rep(" ",n),collapse='')),1,n) }
    lisc2 <- rbind( c("PACKAGE","VERSION","LICENSE"), lisc2)
    lisc2 <- cbind(fixstr(lisc2[,1],36),fixstr(lisc2[,2],15),fixstr(lisc2[,3],30))
    lisc.text <- apply(lisc2, 1, function(s) paste0(s,collapse=''))
    write(lisc.text, "RPackageLicenses.txt")
##  lisc2[grep("LGPL|AGPL|GPL-3",lisc2[,"License"]),]
}
