##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


dataview_module_geneinfo_ui <- function(id, label='', height=c(600,800)) {

    ns <- shiny::NS(id)


    a_OMIM="<a href='https://www.ncbi.nlm.nih.gov/omim/'> OMIM</a>"
    a_KEGG="<a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102409/'> KEGG</a>"
    a_GO="<a href='http://geneontology.org/'>Gene Ontology</a>"

    info_text = paste0('For more information about the the selected gene, follow the hyperlinks to public databases, including ', a_OMIM,', ', a_KEGG, ' and ',a_GO,'.')
        
    PlotModuleUI(
        ns("mod"),
        title = "Gene information",
        label = label,
        outputFunc = htmlOutput,
        outputFunc2 = htmlOutput,        
        info.text = info_text,
        caption = NULL,
        caption2 = NULL,        
        options = NULL,
        download.fmt = NULL,
        width = c("auto","1200"),
        height = height
    )
    
}

dataview_module_geneinfo_server <- function(id,
                                            r.gene = reactive(""),
                                            watermark=FALSE)
{
    moduleServer( id, function(input, output, session) {

        dbg("[dataview_geneinfo_server] created!")
        
        info_data <- shiny::reactive({

            require(org.Hs.eg.db)
            gene <- r.gene()
            req(gene)
            
            gene = toupper(sub(".*:","",gene))
            eg = "1017"
            eg = names(which(as.list(org.Hs.egSYMBOL)==gene))
            eg <- mget(gene, envir=org.Hs.egSYMBOL2EG, ifnotfound=NA)[[1]]
            if(is.na(eg)) eg <- mget(gene, envir=org.Hs.egALIAS2EG, ifnotfound=NA)[[1]]
            eg
            eg = eg[1]
            if(is.null(eg) || length(eg)==0) return(NULL)
            
            res = "(gene info not available)"
            if(length(eg)>0 && !is.na(eg)) {
                info <- getHSGeneInfo(eg)  ## defined in pgx-functions.R
                info$summary <- '(no info available)'
                if(gene %in% names(GENE.SUMMARY)) {
                    info$summary <- GENE.SUMMARY[gene]
                    info$summary <- gsub('Publication Note.*|##.*','',info$summary)
                }

                ## reorder
                nn <- intersect(c("symbol","name","map_location","summary",names(info)),names(info))
                info <- info[nn]
                info$symbol <- paste0(info$symbol,'<br>')
                
                res <- c()
                for(i in 1:length(info)) {
                    xx <- paste(info[[i]], collapse=", ")
                    res[[i]] <- paste0("<b>",names(info)[i],"</b>: ",xx)
                }
                res <- paste(res, collapse="<p>")
            }    
            res
        })


        info.RENDER <- function() {
            res <- info_data()
            shiny::wellPanel(shiny::HTML(res))
        }
        
        modal_info.RENDER <- function() {
            info.RENDER()
        }
        
        PlotModuleServer(
            "mod",
            plotlib = "generic",
            plotlib2 = "generic",
            func = info.RENDER,
            func2 = modal_info.RENDER,
            ##csvFunc = info_data,   ##  *** downloadable data as CSV
            renderFunc = shiny::renderUI,
            renderFunc2 = shiny::renderUI        
        )

    })  ## end of moduleServer
}


