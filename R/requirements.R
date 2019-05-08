##source("http://bioconductor.org/biocLite.R")

if(!require(devtools)) install.packages("devtools")
if(!require(BiocManager)) install.packages("BiocManager")
install.pkg <- function(pkg) {
    if(!pkg %in% installed.packages()) {
        require(BiocManager)
        BiocManager::install(pkg, dependencies=TRUE,
                             ask=FALSE, update=FALSE)
    }
}
install.pkgs <- function(pkgs) {
    for(pkg in pkgs) install.pkg(pkg)
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
##install.pkg("fpc")
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
install.pkg("gridGraphics")
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
install.pkg("NNLM")
install.pkg("nnls")
##install.pkg("clusterProfiler", version = "3.8")
install.pkg("HiveR")

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



