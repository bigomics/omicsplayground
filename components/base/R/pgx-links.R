##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

tagsub <- function(s) {
    s <- gsub("\\{\\{","<i>",gsub("\\}\\}","</i>",s))  ## {{...}}
    s <- gsub("\\{","<code>",gsub("\\}","</code>",s))  ## {...}
    return(s)
}
## example: tagsub("Hi this is an {{option}} block with a {setting}")

a_OMIM="<a href='https://www.ncbi.nlm.nih.gov/omim/'> OMIM</a>"
a_ImmProt="<a href='https://www.ncbi.nlm.nih.gov/pubmed/28263321'> ImmProt</a>"
a_HPA="<a href='https://www.nature.com/articles/nbt1210-1248'> HPA</a>"
a_GTEx="<a href='https://www.ncbi.nlm.nih.gov/pubmed/23715323'> GTEx</a>"

a_PCA="<a href='https://www.ncbi.nlm.nih.gov/pubmed/19377034'> PCA</a>"
a_tSNE="<a href='http://jmlr.org/papers/volume15/vandermaaten14a/vandermaaten14a.pdf'> tSNE</a>"
a_GSVA="<a href='https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7'> GSVA</a>"
a_ssGSEA="<a href='https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7'> ssGSEA</a>"
a_KEGG="<a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102409/'> KEGG</a>"
a_GO="<a href='http://geneontology.org/'>Gene Ontology</a>"
a_MSigDB="<a href='http://software.broadinstitute.org/gsea/msigdb'> GO</a>"
a_Hallmark="<a href='http://software.broadinstitute.org/gsea/msigdb'>Hallmark</a>"
a_Spearman="<a href='https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient'>Spearman rank correlation</a>" 
a_Fisher="<a href='https://www.jstor.org/stable/2340521?seq=1#metadata_info_tab_contents'>Fisher exact test</a>"
a_GSEA="<a href='http://software.broadinstitute.org/gsea/index.jsp'>GSEA</a>" 
a_camera="<a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3458527'>camera</a>"
a_fry="<a href='https://academic.oup.com/bioinformatics/article/26/17/2176/200022'>fry</a>"


SLOGAN = c(
    "'I Feel Empowered' - with OmicsPlayground",
    "'I Love Data' - with OmicsPlayground",
    "'I Get Insights' - with OmicsPlayground",
    "'Data, Knowledge, Insight' - with OmicsPlayground",
    "'So Much More' - with OmicsPlayground",
    "'Do-it-myself' - with OmicsPlayground",
    "'Yes-I-Can' - with OmicsPlayground",
    "'Dig Deeper' - with OmicsPlayground",
    "'Never Stop Exploring' - with OmicsPlayground",
    "'Take control' - with OmicsPlayground",    
    "'My Eureka! moments' - with OmicsPlayground",
    "'I See Clearly Now' - with OmicsPlayground"
    )
