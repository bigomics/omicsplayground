##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


dataview_module_geneinfo_ui <- function(
    id,
    label = "",
    title,
    height,
    width,
    caption,
    info.text,
    info.methods,
    info.references) {
  ns <- shiny::NS(id)


  a_OMIM <- "<a href='https://www.ncbi.nlm.nih.gov/omim/'> OMIM</a>"
  a_KEGG <- "<a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102409/'> KEGG</a>"
  a_GO <- "<a href='http://geneontology.org/'>Gene Ontology</a>"

  PlotModuleUI(
    ns("mod"),
    title = title,
    label = label,
    outputFunc = htmlOutput,
    outputFunc2 = htmlOutput,
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    caption = caption,
    caption2 = NULL,
    options = NULL,
    download.fmt = NULL,
    width = width,
    height = height
  )
}

dataview_module_geneinfo_server <- function(id,
                                            pgx,
                                            r.gene = reactive(""),
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    geneinfo_data <- shiny::reactive({
      gene <- r.gene()
      shiny::req(gene)

      jj <- match(gene, rownames(pgx$genes))
      symbol <- pgx$genes$symbol[jj]
      info <- playbase::getHSGeneInfo(symbol) ## defined in pgx-functions.R
      res <- tspan("(gene info not available)")
      if (!is.null(info)) {
        info$summary <- "(no info available)"
        if (symbol %in% names(playdata::GENE_SUMMARY)) {
          info$summary <- playdata::GENE_SUMMARY[symbol]
          info$summary <- gsub("Publication Note.*|##.*", "", info$summary)
        }

        ## reorder
        nn <- intersect(
          c("symbol", "name", "map_location", "summary", names(info)),
          names(info)
        )
        info <- info[nn]
        info$symbol <- paste0(info$symbol, "<br>")

        res <- c()
        for (i in 1:length(info)) {
          xx <- paste(info[[i]], collapse = ", ")
          res[[i]] <- paste0("<b>", names(info)[i], "</b>: ", xx)
        }
        res <- paste(res, collapse = "<p>")
      }
      res
    })


    info.RENDER <- function() {
      res <- geneinfo_data()
      div(shiny::HTML(res), class = "gene-info")
    }

    modal_info.RENDER <- function() {
      res <- geneinfo_data()
      div(shiny::HTML(res), class = "gene-info", style = "font-size:1.3em;")
    }

    PlotModuleServer(
      "mod",
      plotlib = "generic",
      plotlib2 = "generic",
      func = info.RENDER,
      func2 = modal_info.RENDER,
      ## csvFunc = geneinfo_data,   ##  *** downloadable data as CSV
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI
    )
  }) ## end of moduleServer
}

## dataview_module_geneinfo_server <- function(id,
##                                             r.gene = reactive(""),
##                                             watermark = FALSE) {
##   moduleServer(id, function(input, output, session) {
##     geneinfo_data <- shiny::reactive({
##       gene <- r.gene()
##       req(gene)

##       gene <- toupper(sub(".*:", "", gene))
##       eg <- AnnotationDbi::mget(gene,
##         envir = org.Hs.eg.db::org.Hs.egSYMBOL2EG,
##         ifnotfound = NA
##       )[[1]]
##       if (isTRUE(is.na(eg))) {
##         eg <- AnnotationDbi::mget(gene,
##           envir = org.Hs.eg.db::org.Hs.egALIAS2EG, ifnotfound = NA
##         )[[1]]
##       }
##       eg <- eg[1]
##       if (is.null(eg) || length(eg) == 0) {
##         return(NULL)
##       }

##       res <- "(gene info not available)"
##       if (length(eg) > 0 && !is.na(eg)) {
##         info <- playbase::getHSGeneInfo(eg) ## defined in pgx-functions.R
##         info$summary <- "(no info available)"
##         if (gene %in% names(playdata::GENE_SUMMARY)) {
##           info$summary <- playdata::GENE_SUMMARY[gene]
##           info$summary <- gsub("Publication Note.*|##.*", "", info$summary)
##         }

##         ## reorder
##         nn <- intersect(c("symbol", "name", "map_location", "summary", names(info)), names(info))
##         info <- info[nn]
##         info$symbol <- paste0(info$symbol, "<br>")

##         res <- c()
##         for (i in 1:length(info)) {
##           xx <- paste(info[[i]], collapse = ", ")
##           res[[i]] <- paste0("<b>", names(info)[i], "</b>: ", xx)
##         }
##         res <- paste(res, collapse = "<p>")
##       }
##       res
##     })


##     info.RENDER <- function() {
##       res <- geneinfo_data()
##       div(shiny::HTML(res), class = "gene-info")
##     }

##     modal_info.RENDER <- function() {
##       res <- geneinfo_data()
##       div(shiny::HTML(res), class = "gene-info", style = "font-size:1.3em;")
##     }

##     PlotModuleServer(
##       "mod",
##       plotlib = "generic",
##       plotlib2 = "generic",
##       func = info.RENDER,
##       func2 = modal_info.RENDER,
##       ## csvFunc = geneinfo_data,   ##  *** downloadable data as CSV
##       renderFunc = shiny::renderUI,
##       renderFunc2 = shiny::renderUI
##     )
##   }) ## end of moduleServer
## }
