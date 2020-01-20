library(shiny)
library(shinyjs)
library(shinyBS)
library(shinyjqui)
library(rmarkdown)
library(shinycssloaders)
library(dragulaR)
library(shinyWidgets)
library(pryr)

library(survival)
library(knitr)
library(scatterD3)
library(fastcluster)
library(ComplexHeatmap)
library(plotly)
library(Matrix)
library(igraph)
library(DT)
library(ggplot2)
library(data.table)
library(dplyr)
library(org.Hs.eg.db)
## library(Cairo)
library(dragulaR)

    
##useShinyjs(rmd=TRUE)  
useShinyjs()
ht_global_opt(fast_hclust = TRUE)
options(shiny.maxRequestSize = 500*1024^2)  ## max 200Mb upload

dbg <- function(... ) {
    if(DEV.VERSION) {
        ##msg = paste0(ifelse(is.null(module),"",paste0("<",module,"> ")),msg)
        msg = sapply( list(...),paste,collapse=" ")
        cat(paste0("DBG ",sub("\n$","",paste(msg,collapse=" ")),"\n"))
    }
}

source(file.path(RDIR,"gx-heatmap.r"))
source(file.path(RDIR,"gx-plot.r"))
source(file.path(RDIR,"gx-limma.r"))
source(file.path(RDIR,"gx-volcano.r"))
source(file.path(RDIR,"gx-combat.r"))
source(file.path(RDIR,"gx-util.r"))
       
source(file.path(RDIR,"gset-gsea.r"))
source(file.path(RDIR,"gset-fisher.r"))
source(file.path(RDIR,"gset-meta.r"))

source(file.path(RDIR,"ngs-cook.r"))
source(file.path(RDIR,"ngs-fit.r"))
source(file.path(RDIR,"ngs-functions.R"))

source(file.path(RDIR,"pgx-functions.R"))
source(file.path(RDIR,"pgx-contrasts.R"))
source(file.path(RDIR,"pgx-graph.R"))
source(file.path(RDIR,"pgx-deconv.R"))
source(file.path(RDIR,"pgx-cna.R"))
source(file.path(RDIR,"pgx-plotting.R"))
source(file.path(RDIR,"pgx-correct.R"))
source(file.path(RDIR,"pgx-predict.R"))
source(file.path(RDIR,"pgx-links.R"))
source(file.path(RDIR,"pgx-modules.R"))
source(file.path(RDIR,"pgx-upload.R"))
source(file.path(RDIR,"pgx-proteomics.R"))
source(file.path(RDIR,"pgx-drugs.R"))
source(file.path(RDIR,"pgx-files.R"))
source(file.path(RDIR,"pgx-wordcloud.R"))
source(file.path(RDIR,"pgx-ui.R"))
source(file.path(RDIR,"pgx-tcga.R"))
source(file.path(RDIR,"pgx-archs4.R"))
source(file.path(RDIR,"pgx-cluster.R"))
source(file.path(RDIR,"pgx-signature.R"))

source(file.path(RDIR,"xcr-graph.r"))
source(file.path(RDIR,"ui-code.R"))

source(file.path(RDIR,"compute2-genes.R"))
source(file.path(RDIR,"compute2-genesets.R"))
source(file.path(RDIR,"compute2-extra.R"))

