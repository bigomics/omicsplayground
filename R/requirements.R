##
##
## NOTE: This file is supposed to run in the folder .../R/
##

if(!require(devtools)) install.packages("devtools")
install.packages("BiocManager", version="3.10")
require(devtools)
require(BiocManager)

LOCAL.PKGS <- sub("_.*","",dir("../ext/packages"))
LOCAL.PKGS

install.pkg <- function(pkg, force=FALSE) {
    if(force || !pkg %in% installed.packages()) {
        if(pkg %in% LOCAL.PKGS) {
            ## if available locally, we install local version
            cat("installing",pkg,"from local folder...\n")
            pkg1 = dir("../ext/packages",pattern=paste0(pkg,"_"),full.name=TRUE)
            try(install.packages(pkg1,repos=NULL,type="source"))
        } else {
            cat("installing",pkg,"from CRAN/BioConductor...\n")
            try(BiocManager::install(pkg, dependencies=NA,
                                      ask=FALSE, update=FALSE))
        }
    } else {
        cat("package",pkg,"already installed\n")
    }
}
install.pkgs <- function(pkgs, force=FALSE) {
    for(pkg in pkgs) install.pkg(pkg, force=force)
}
remove.pkg <- function(pkg) {
    if(pkg %in% installed.packages()) remove.packages(pkg)
}
remove.pkgs <- function(pkgs, force=FALSE) {
    for(pkg in pkgs) remove.pkg(pkg)
}

BIG.NOTUSED <- c(
    "reactome.db", ## >2GB!!
    "BH","PCSF",
    "DeMixT", ## purify
    "RNAseqData.HNRNPC.bam.chr14",
    ##"org.Mm.eg.db",
    "tximportData",
    ##"EnsDb.Hsapiens.v86",
    ##"EnsDb.Mmusculus.v79",
    "TxDb.Hsapiens.UCSC.hg19.knownGene",
    "TxDb.Mmusculus.UCSC.mm10.knownGene"
)    

PKG.MANUAL <- c(
    "gputools","Seurat","EPIC","PCSF","NNLM","iTALK",
    "fpc","grid","gridGraphics","Rgraphviz","rWordCloud",
    "shinyparticles","FastGGM","monocle3","proteus",
    ## "infercnv","pathview",
    "mygene","diptest")

##---------------------------------------------------------------------
## Install base packages
##---------------------------------------------------------------------

base.pkg = c("shiny","flexdashboard","shinydashboard",
             "shinydashboardPlus",'R.utils')
install.pkgs(base.pkg)

##---------------------------------------------------------------------
## Automatically scan all used packages and install
##---------------------------------------------------------------------

pkg.used <- system("grep 'library(\\|require(' *R *r ../shiny/*R ../shiny/modules/*R", intern=TRUE)
pkg.used <- gsub(".*require\\(|.*library\\(","",pkg.used)
pkg.used <- gsub("\"|\'|\\).*","",pkg.used)
pkg.used <- grep("[ ]|quietly",pkg.used,value=TRUE,invert=TRUE)

pkg.needed <- c('umap','corrplot','wordcloud',"optparse","docopt","randomForest")
pkg.used <- c(pkg.used, pkg.needed)
pkg.used <- sort(unique(pkg.used))

install.pkgs( setdiff(pkg.used, c(PKG.MANUAL,BIG.NOTUSED)) )

##---------------------------------------------------------------------
## reinstall problematics ones
##---------------------------------------------------------------------

require(devtools)
##install.pkg("grid", force=TRUE)
install.pkgs(c("gridGraphics","Rgraphviz","fastcluster", "mygene",
               "diptest", "fpc", "webshot"))
##webshot::install_phantomjs() 
devtools::install_version("mnormt", version="1.5-7", repos="http://cran.us.r-project.org") 
install.pkgs(c('umap','corrplot','wordcloud','metap','brew','Seurat'))

##---------------------------------------------------------------------
## Install latest from GITHUB
##---------------------------------------------------------------------
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
##devtools::install_github("IOR-Bioinformatics/PCSF", dependencies=TRUE, type="source")
devtools::install_github('linxihui/NNLM')
devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)
devtools::install_github('adymimos/rWordCloud', force=TRUE)
## remotes::install_github("dreamRs/shinyparticles")
remotes::install_github("dreamRs/particlesjs")
remotes::install_github("JohnCoene/waiter")

##---------------------------------------------------------------------
## ONLY DEV.VERSION (single-cell trajectories)
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

install.pkgs(c("RcppParallel"))
devtools::install_github('wt2015-github/FastGGM')
install.pkgs(c("HiddenMarkov","coin","rjags","future","argparse"))
install.pkg("rjags", force=TRUE)
##install.packages("../ext/packages/infercnv_1.1.3mod.tar.gz",repos=NULL,type="source")  ## old version
BiocManager::install("infercnv")
install.pkgs(c("KEGGREST","pathview"), force=TRUE)

##---------------------------------------------------------------------
## remove unneccessary big packages...
##---------------------------------------------------------------------
BIG.NOTUSED
remove.pkgs(BIG.NOTUSED)
