##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##























## 

    
##useShinyjs(rmd=TRUE)  
shinyjs::useShinyjs()
ComplexHeatmap::ht_global_opt(fast_hclust = TRUE)
options(shiny.maxRequestSize = 500*1024^2)  ## max 200Mb upload

dbg <- function(... ) {
    ##msg = paste0(ifelse(is.null(module),"",paste0("<",module,"> ")),msg)
    msg = sapply( list(...),paste,collapse=" ")
    message(paste0("DBG ",sub("\n$","",paste(msg,collapse=" "))))
}

##source(file.path(RDIR,"pgx-functions.R"), local=TRUE)
source(file.path(RDIR,"pgx-functions.R"))
##source(file.path(RDIR,"pgx-init.R"))

source(file.path(RDIR,"gx-heatmap.r"), local=TRUE)
source(file.path(RDIR,"gx-plot.r"))
source(file.path(RDIR,"gx-limma.r"))
source(file.path(RDIR,"gx-volcano.r"))
source(file.path(RDIR,"gx-combat.r"))
source(file.path(RDIR,"gx-util.r"))
       
source(file.path(RDIR,"gset-gsea.r"))
source(file.path(RDIR,"gset-fisher.r"))
source(file.path(RDIR,"gset-rankcor.r"))
source(file.path(RDIR,"gset-meta.r"))

source(file.path(RDIR,"ngs-cook.r"))
source(file.path(RDIR,"ngs-fit.r"))

source(file.path(RDIR,"pgx-api.R"), local=TRUE)
source(file.path(RDIR,"pgx-contrasts.R"))
source(file.path(RDIR,"pgx-graph.R"))
source(file.path(RDIR,"pgx-deconv.R"))
source(file.path(RDIR,"pgx-cna.R"))
source(file.path(RDIR,"pgx-plotting.R"))
source(file.path(RDIR,"pgx-correct.R"))
source(file.path(RDIR,"pgx-predict.R"))
source(file.path(RDIR,"pgx-links.R"))
source(file.path(RDIR,"pgx-modules.R"), local=TRUE)
##source(file.path(RDIR,"pgx-modules.R"))
source(file.path(RDIR,"pgx-proteomics.R"))
source(file.path(RDIR,"pgx-drugs.R"))
source(file.path(RDIR,"pgx-files.R"))
source(file.path(RDIR,"pgx-wordcloud.R"))
source(file.path(RDIR,"pgx-ui.R"))
source(file.path(RDIR,"pgx-tcga.R"))
source(file.path(RDIR,"pgx-cluster.R"))
source(file.path(RDIR,"pgx-signature.R"))
##source(file.path(RDIR,"pgx-archs4.R"))
##source(file.path(RDIR,"pgx-getgeo.R"))
source(file.path(RDIR,"pgx-correlation.R"))

source(file.path(RDIR,"xcr-graph.r"))
source(file.path(RDIR,"ui-code.R"))

source(file.path(RDIR,"compute2-genes.R"))
source(file.path(RDIR,"compute2-genesets.R"))
source(file.path(RDIR,"compute2-extra.R"))
source(file.path(RDIR,"pgx-compute.R"), local=TRUE)
##source(file.path(RDIR,"pgx-compute.R"))
source(file.path(RDIR,"pgx-singlecell.R"))
source(file.path(RDIR,"pgx-vizpanels.R"))


