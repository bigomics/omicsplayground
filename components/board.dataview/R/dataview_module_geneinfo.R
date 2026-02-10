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
  info.references
) {
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

      feature <- r.gene()
      shiny::req(feature %in% rownames(pgx$genes))

      info <- playbase::pgx.getFeatureInfo(pgx, feature)

      if (is.null(info)) {
        info <- tspan("(gene info not available)")
        return(info)
      }

      names(info) <- tolower(names(info))
      names(info) <- sub("gene_symbol", "symbol", names(info))
      names(info) <- sub("gene_title", "title", names(info))
      names(info) <- sub("uniprot", "protein", names(info))
      names(info) <- sub("map_location", "genome location", names(info))
      names(info) <- sub("databases", "links", names(info))

      nn1 <- intersect(
        c(
          "feature", "gene_symbol", "symbol",
          "name", "gene_name",
          "organism", "human_ortholog", "ortholog",
          "datatype", "data_type",
          "uniprot", "protein",
          "map_location", "map", "genome location",
          "databases", "links"
        ),
        names(info)
      )
      nn2 <- intersect(
        c(
          "title", "gene_title",
          "summary", "description"
        ),
        names(info)
      )
      nn3 <- setdiff(names(info), c(nn1, nn2))
      info <- info[c(nn1, nn2, nn3)]

      res <- c()
      for (i in 1:length(info)) {
        xx <- paste(info[[i]], collapse = ", ")
        res[[i]] <- paste0("<b>", names(info)[i], "</b>: ", xx)
      }
      names(res) <- names(info)
      res <- c(
        "<p>", paste(res[nn1], collapse = "<br>"),
        "<p>", paste(res[nn2], collapse = "<p>"),
        "<p>", paste(res[nn3], collapse = "<br>")
      )

      return(res)

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
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI
    )
  })
}
