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
    ## prepare data
    geneinfo_data <- shiny::reactive({
      feature <- r.gene()
      shiny::req(feature %in% rownames(pgx$X))

      organism <- pgx$organism
      if (is.null(organism) || is.na(organism) || organism == "") organism <- "Human"
      datatype <- pgx$datatype
      if (is.null(datatype) || is.na(datatype) || datatype == "") datatype <- "RNA-seq"

      jj <- match(feature, rownames(pgx$genes))
      symbol <- pgx$genes$symbol[jj]
      ortholog <- pgx$genes$human_ortholog[jj]

      if (datatype == "metabolomics") {
        info <- playbase::getMetaboliteInfo(
          organism = organism,
          chebi = symbol
        )
      } else {
        info <- playbase::getOrgGeneInfo(
          organism = organism,
          gene = symbol,
          feature = feature,
          ortholog = ortholog,
          datatype = datatype,
          as.link = TRUE
        )

        ## reorder
        nn <- intersect(
          c(
            "gene_symbol", "organism", "name", "map_location",
            "uniprot", "databases", "summary", names(info)
          ),
          names(info)
        )

        info <- info[nn]
        names(info) <- sub("gene_symbol", "symbol", names(info))
        names(info) <- sub("uniprot", "protein", names(info))
        names(info) <- sub("map_location", "genome location", names(info))
      }

      if (is.null(info)) {
        info$summary <- "(no info available)"
        if (symbol %in% names(playdata::GENE_SUMMARY)) {
          info$summary <- playdata::GENE_SUMMARY[symbol]
          info$summary <- gsub("Publication Note.*|##.*", "", info$summary)
        }
      }

      ## add feature name is not symbol
      info$feature <- NULL
      if (!is.na(symbol) && feature != symbol) {
        info <- c(feature = feature, info)
      }
      ## info$organism <- NULL
      ## info$databases <- NULL

      # prepare info for display
      res <- c()
      for (i in 1:length(info)) {
        xx <- paste(info[[i]], collapse = ", ")
        res[[i]] <- paste0("<b>", names(info)[i], "</b>: ", xx)
      }
      res <- paste(res, collapse = "<p>")

      if (is.null(res)) {
        res <- tspan("(gene info not available)")
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
