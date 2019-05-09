##source("http://bioconductor.org/biocLite.R")

install.packages("devtools")
install.packages("BiocManager")
require(devtools)
require(BiocManager)

install.pkg <- function(pkg, force=FALSE) {
    if(force && (pkg %in% installed.packages())) {
        ##remove.packages(pkg)
    }
    if(force || !pkg %in% installed.packages()) {
        require(BiocManager)
        cat("installing",pkg,"\n")
        BiocManager::install(pkg, dependencies=NA,
                             ask=FALSE, update=FALSE)
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


##---------------------------------------------------------------------
## CRAN packages
##---------------------------------------------------------------------

install.pkg("Rcpp")
install.pkg("XML")
install.pkg("shinyjs")
install.pkg("flexdashboard")
install.pkg("shinyWidgets")
install.pkg("kableExtra")

install.pkg("DT")
install.pkg("htmltools")
install.pkg("reticulate")
install.pkg("tidyverse")
install.pkg("qlcMatrix")
install.pkg("ComplexHeatmap")
install.pkg("gmodels")
install.pkg("matrixTests")
install.pkg("metap")
install.pkg("plotly")
install.pkg("dplyr")
install.pkg("scatterD3")
install.pkg("threejs")
install.pkg("locfit")
install.pkg("irlba")
install.pkg("ggpubr")
install.pkg("corrplot")
install.pkg("corpora")  ## for fisher.pval
install.pkg("xgboost")
install.pkg("randomForest")
install.pkg("randomForestSRC")
install.pkg("rpart")
install.pkg("rpart.plot")
install.pkg("party")
install.pkg("partykit")
##install.pkg("NNLM", force=TRUE)
install_version("NNLM", version="0.4.1", repos="http://cran.us.r-project.org")
install.pkg("nnls")
install.pkg("HiveR")
install.pkg("grid", force=TRUE)

## problematics ones
require(devtools)
##install.pkg("fpc", force=TRUE)
##install.packages("gridGraphics")
install_version("gridGraphics", version="0.3-0", repos="http://cran.us.r-project.org")
install.pkg("fastcluster", force=TRUE)

##---------------------------------------------------------------------
## Bioconductor packages
##---------------------------------------------------------------------
install.pkg("preprocessCore")
install.pkg("BiocParallel")
install.pkg("BiocGenerics")
install.pkg("org.Hs.eg.db")
install.pkg("org.Mm.eg.db")
install.pkg("EnsDb.Hsapiens.v86")
##install.pkg("EnsDb.Mmusculus.v79")
install.pkg("hgu133plus2.db")

install.pkg("GSVA")
install.pkg("fgsea")
install.pkg("sva")
install.pkg("tximport")
install.pkg("limma")
install.pkg("edgeR")
install.pkg("DESeq2")
install.pkg("ensembldb")
install.pkg("pcaMethods")
install.pkg("DeconRNASeq")
install.pkg("GenomicRanges")
install.pkg("IRanges")
install.pkg("KEGG.db")
install.pkg("KEGGREST")
install.pkg("KEGGgraph")
install.pkg("GO.db")
##install.pkg("pathview")
install.pkg("ComICS")
install.pkg("scran")
install.pkg("SingleCellExperiment")
install.pkg("SummarizedExperiment")
install.pkg("diffusionMap")
install.pkg("topGO")
install.pkg("mixOmics")
install.pkg("mygene")
install.pkg("Rgraphviz", force=TRUE)


