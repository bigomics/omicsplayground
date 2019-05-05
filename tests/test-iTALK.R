## This example data is from 10x pbmc dataset. Samples are randomly
## selected from each cell type. And groups are randomly assigned to
## each sample to make the illustration.

library(iTALK)

# read the data
data <- read.table('~/Downloads/iTALK/example/example_data.txt',sep='\t',header=T,stringsAsFactors = F)
highly_exprs_genes <- rawParse(data,top_genes=50,stats='mean')

## find the ligand-receptor pairs from highly expressed genes
comm_type='growth factor'
comm_type='checkpoint'
comm_list <- c('growth factor','other','cytokine','checkpoint')
cell_col <- structure(c('#4a84ad','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b','#9ef49a'),names=unique(data$cell_type))
par(mfrow=c(1,2))
res<-NULL
for(comm_type in comm_list){

    res_cat <- FindLR(highly_exprs_genes,datatype='mean count',comm_type=comm_type)
    res_cat <- res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
    #plot by ligand category
    #overall network plot
    NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)

    ##top 20 ligand-receptor pairs
    res_cat[1:20,]
    LRPlot(res_cat[1:20,], datatype='mean count', cell_col=cell_col,
           link.arr.lwd=res_cat$cell_from_mean_exprs[1:20],
           link.arr.width=res_cat$cell_to_mean_exprs[1:20])
    title(comm_type)
    res<-rbind(res,res_cat)

}
res<-res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:20,]
NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
LRPlot(res[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res$cell_from_mean_exprs[1:20],link.arr.width=res$cell_to_mean_exprs[1:20])


##-------------------------------------------------------------------------------
## significant ligand-receptor pairs between compare groups
##-------------------------------------------------------------------------------
require(dplyr)

# randomly assign the compare group to each sample
data <- data %>% mutate(compare_group=sample(2,nrow(data),replace=TRUE))

## find DEGenes of regulatory T cells and NK cells between these 2 groups
deg_t  <- DEG(data %>% filter(cell_type=='regulatory_t'),method='Wilcox',contrast=c(2,1))
deg_nk <- DEG(data %>% filter(cell_type=='cd56_nk'),method='Wilcox',contrast=c(2,1))

## find significant ligand-receptor pairs and do the plotting
par(mfrow=c(1,2))
res<-NULL
for(comm_type in comm_list){
    res_cat<-FindLR(deg_t,deg_nk,datatype='DEG',comm_type=comm_type)
    res_cat<-res_cat[order(res_cat$cell_from_logFC*res_cat$cell_to_logFC,decreasing=T),]
    ##plot by ligand category
    if(nrow(res_cat)==0){
        next
    } else if(nrow(res_cat>=20)){
        LRPlot(res_cat[1:20,],datatype='DEG',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_logFC[1:20],link.arr.width=res_cat$cell_to_logFC[1:20])
    } else {
        LRPlot(res_cat,datatype='DEG',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_logFC,link.arr.width=res_cat$cell_to_logFC)
    }
    NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
    title(comm_type)
    res<-rbind(res,res_cat)
}
if(is.null(res)){
    print('No significant pairs found')
} else if(nrow(res)>=20){
    res<-res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),][1:20,]
    NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
    LRPlot(res[1:20,],datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC[1:20],link.arr.width=res$cell_to_logFC[1:20])
} else {
    NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
    LRPlot(res,datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC,link.arr.width=res$cell_to_logFC)
}




BiocManager::install("scTensor")
require(scTensor)
